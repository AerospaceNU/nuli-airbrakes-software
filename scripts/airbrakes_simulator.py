from rocketpy import Environment, Rocket, SolidMotor, Flight

"""
No motor added in simulation as this is simulating post-burnout
"""


class AirbrakesSimulator:
    def __init__(self, launch_altitude=0, latitude=44.826255, longitude=-73.177421):
        self.env = Environment(
            latitude=latitude,
            longitude=longitude,
            elevation=launch_altitude
        )

        # rocket components
        self.rocket = None
        self.nose = None
        self.fins = None
        self.tail = None
        self.parachute_main = None
        self.parachute_drogue = None
        self.airbrakes = None

        # set atmosphere, or use real weather data
        self.env.set_atmospheric_model(type="standard_atmosphere")

    def create_rocket(self,
                      radius,  # m
                      length,  # m (total rocket length)
                      nose_cone_length,  # m
                      dry_mass,  # kg (without motor)
                      inertia,  # (Ixx, Iyy, Izz) in kg*m^2
                      center_of_mass_without_motor,
                      coordinate_system_orientation="tail_to_nose",
                      power_off_drag="rocketpy-data/data/rockets/calisto/powerOffDragCurve.csv",  #
                      power_on_drag="rocketpy-data/data/rockets/calisto/powerOnDragCurve.csv"):

        # Estimate moments of inertia if not provided (thin rod approximation)
        # if inertia[0] is None:
        #     # Longitudinal (roll) inertia
        #     i_xx = 0.5 * dry_mass * (radius)**2
        #     # Lateral (pitch/yaw) inertia
        #     i_yy = i_zz = (1/12) * dry_mass * length**2 + dry_mass * (center_of_mass - length/2)**2
        #     inertia = (i_xx, i_yy, i_zz)

        # Rocket body
        self.rocket = Rocket(
            radius=radius,
            mass=dry_mass,
            inertia=inertia,
            power_on_drag=power_on_drag,
            power_off_drag=power_off_drag,
            center_of_mass_without_motor=center_of_mass_without_motor,
            coordinate_system_orientation=coordinate_system_orientation
        )

        # nose cone
        self.rocket.add_nose(length=nose_cone_length,
                             kind="ogive",
                             position=1.278  # because tail to nose
                             )

        # fins
        self.rocket.add_trapezoidal_fins(
            n=4,
            root_chord=0.120,
            tip_chord=0.040,
            span=0.100,
            sweep_length=None,
            cant_angle=0,
            position=length * 0.9,  # near the base
        )

        # tail
        self.rocket.add_tail(
            top_radius=0.0635, bottom_radius=0.0435, length=0.060, position=-1.194656
        )

        # parachute - main
        self.rocket.add_parachute(
            name="main",
            cd_s=10.0,
            trigger=800,  # ejection altitude in meters
            sampling_rate=105,
            lag=1.5,
            noise=(0, 8.3, 0.5),
        )

        # parachute - drogue
        self.rocket.add_parachute(
            name="drogue",
            cd_s=1.0,
            trigger="apogee",  # ejection at apogee
            sampling_rate=105,
            lag=1.5,
            noise=(0, 8.3, 0.5),
        )

    def add_airbrakes(self,
                      airbrake_area,  # m^2 per panel
                      num_panels=4,
                      drag_coefficient_deployed=1.2,
                      deployment_altitude=100,  # m AGL
                      deployment_velocity_threshold=10):  # m/s

        # Calculate additional drag area from airbrakes
        total_airbrake_area = airbrake_area * num_panels

        def controller_function(
                time, sampling_rate, state, state_history, observed_variables, air_brakes
        ):
            # state = [x, y, z, vx, vy, vz, e0, e1, e2, e3, wx, wy, wz]
            altitude_ASL = state[2]
            altitude_AGL = altitude_ASL - self.env.elevation
            vx, vy, vz = state[3], state[4], state[5]

            # Get winds in x and y directions
            wind_x, wind_y = self.env.wind_velocity_x(altitude_ASL), self.env.wind_velocity_y(altitude_ASL)

            # Calculate Mach number
            free_stream_speed = (
                                        (wind_x - vx) ** 2 + (wind_y - vy) ** 2 + (vz) ** 2
                                ) ** 0.5
            mach_number = free_stream_speed / self.env.speed_of_sound(altitude_ASL)

            # Get previous state from state_history
            previous_state = state_history[-1]
            previous_vz = previous_state[5]

            # If we wanted to we could get the returned values from observed_variables:
            # returned_time, deployment_level, drag_coefficient = observed_variables[-1]

            # # Check if the rocket has reached burnout
            # if time < Pro75M1670.burn_out_time:
            #     return None

            # If below 1500 meters above ground level, air_brakes are not deployed
            if altitude_AGL < 1500:
                air_brakes.deployment_level = 0

            # Else calculate the deployment level
            else:
                # Controller logic
                new_deployment_level = (
                        air_brakes.deployment_level + 0.1 * vz + 0.01 * previous_vz ** 2
                )

                # Limiting the speed of the air_brakes to 0.2 per second
                # Since this function is called every 1/sampling_rate seconds
                # the max change in deployment level per call is 0.2/sampling_rate
                max_change = 0.2 / sampling_rate
                lower_bound = air_brakes.deployment_level - max_change
                upper_bound = air_brakes.deployment_level + max_change
                new_deployment_level = min(max(new_deployment_level, lower_bound), upper_bound)

                air_brakes.deployment_level = new_deployment_level

            # Return variables of interest to be saved in the observed_variables list
            return (
                time,
                air_brakes.deployment_level,
                air_brakes.drag_coefficient(air_brakes.deployment_level, mach_number),
            )

        def airbrake_trigger(time, state_vector):
            """
            Trigger function for airbrake deployment
            state_vector contains: x, y, z, vx, vy, vz, e0, e1, e2, e3, omega1, omega2, omega3
            """
            # Get altitude (z position)
            altitude = state_vector[2]

            # Get vertical velocity
            vz = state_vector[5]

            # Deploy if above deployment altitude and still ascending
            if altitude >= deployment_altitude and vz > deployment_velocity_threshold:
                return drag_coefficient_deployed
            else:
                return 0  # Not deployed

        self.airbrakes = self.rocket.add_air_brakes(
            drag_coefficient_curve="./rocketpy-data/data/rockets/calisto/air_brakes_cd.csv",
            controller_function=airbrake_trigger,
            sampling_rate=10,  # Hz
            reference_area=total_airbrake_area,
            clamp=True,
            initial_observed_variables=[0, 0, 0],
            override_rocket_drag=False,
            name="Air Brakes",
        )

    def simulate_flight(
            self,
            launch_angle=90,  # degrees from horizontal
            heading=0,  # degrees from North

    ):
        flight = Flight(
            rocket=self.rocket,
            environment=self.env,
            rail_length=0.01,  # no rail length for post burnout flight
            inclination=launch_angle,
            heading=heading,
            # Airbrakes specific options
            # To simulate the air brakes successfully, we must set
            # time_overshoot to False. This way the simulation will run at the
            # time step defined by our controller sampling rate. Be aware that
            # this will make the simulation run much slower.
            time_overshoot=False,
            terminate_on_apogee=True  # no need for simulation after apogee
        )

        return flight


if __name__ == "__main__":
    simulation = AirbrakesSimulator(launch_altitude=0)

    # create rocket
    simulation.create_rocket(
        radius=0.0635,   # m
        length=1.5,     # m
        nose_cone_length=0.5,   # m
        dry_mass=5.0,   # kg
        center_of_mass_without_motor=0.75,
        inertia=(6.321, 6.321, 0.034)
    )

    simulation.add_airbrakes(
        airbrake_area=1,    # m^2 per panel
        num_panels=4,
        drag_coefficient_deployed=1.28,
        deployment_altitude=100,  # m
    )

    flight_sim = simulation.simulate_flight()
    flight_sim.info()

"""
    dry_mass=5.0,  # kg
    diameter=0.15,  # m
    length=1.5,  # m
    initial_velocity=150,  # m/s
    airbrake_area=0.002,  # m^2 per panel
    deployment_altitude=100  # m
"""