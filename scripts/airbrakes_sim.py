from rocketpy import Environment, SolidMotor, Rocket, Flight
import matplotlib.pyplot as plt


class AirbrakesSimulation:
    def __init__(self):
        self.env = Environment(latitude=32.990254, longitude=-106.974998, elevation=1400)

        self.AeroTech_M2100G = SolidMotor(
            thrust_source="./rocketpy-data/data/motors/aerotech/AeroTech_M2100G.eng",
            dry_mass=2.503,  # kg -- 5
            dry_inertia=(0.125, 0.125, 0.002),
            nozzle_radius=33 / 1000,
            grain_number=5,
            grain_density=1815,
            grain_outer_radius=33 / 1000,
            grain_initial_inner_radius=15 / 1000,
            grain_initial_height=120 / 1000,
            grain_separation=5 / 1000,
            grains_center_of_mass_position=0.397,
            center_of_dry_mass_position=0.317,
            nozzle_position=0,
            burn_time=3.541,
            throat_radius=11 / 1000,
            coordinate_system_orientation="nozzle_to_combustion_chamber",
        )

        self.calisto = Rocket(
            radius=127 / 2000,
            mass=11.732,
            inertia=(6.13, 7.13, 0.034),  # (6.321, 6.321, 0.034)
            power_off_drag="./rocketpy-data/data/rockets/calisto/powerOffDragCurve.csv",
            power_on_drag="./rocketpy-data/data/rockets/calisto/powerOnDragCurve.csv",
            center_of_mass_without_motor=1.38,
            coordinate_system_orientation="nose_to_tail",
        )

        self.rail_buttons = self.calisto.set_rail_buttons(
            upper_button_position=0.0818,
            lower_button_position=-0.618,
            angular_position=45,
        )

        self.calisto.add_motor(self.AeroTech_M2100G, position=-1.255)

        self.nose_cone = self.calisto.add_nose(
            length=0.686, kind="ogive", position=0.686
        )

        self.fin_set = self.calisto.add_trapezoidal_fins(
            n=3,
            root_chord=0.152,
            tip_chord=0.102,
            span=0.127,
            position=2.295,
            cant_angle=0.5,
            airfoil=("./rocketpy-data/data/airfoils/NACA0012-radians.txt", "radians"),
        )

        self.tail = self.calisto.add_tail(
            top_radius=0.0635, bottom_radius=0.0636, length=0.914, position=1.55
        )

    #### AIRBRAKES
    def airbrakes_controller_function(
            self, time, sampling_rate, state, state_history, observed_variables, air_brakes
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

        # Check if the rocket has reached burnout
        if time < self.AeroTech_M2100G.burn_out_time:
            return None

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

    def add_airbrakes(
            self,
            num_airbrakes,  # int
            airbrakes_area  # m^2
    ):
        airbrakes_reference_area = num_airbrakes * airbrakes_area

        air_brakes = self.calisto.add_air_brakes(
            drag_coefficient_curve="./rocketpy-data/data/rockets/calisto/air_brakes_cd.csv",
            controller_function=self.airbrakes_controller_function,
            sampling_rate=10,
            reference_area=airbrakes_reference_area,
            clamp=True,
            initial_observed_variables=[0, 0, 0],
            override_rocket_drag=False,
            name="Air Brakes",
        )

        return air_brakes


if __name__ == "__main__":
    # Normal w/o Airbrakes
    sim_normal = AirbrakesSimulation()
    test_flight_normal = Flight(
        rocket=sim_normal.calisto,
        environment=sim_normal.env,
        rail_length=5.2,
        inclination=85,
        heading=0,
        time_overshoot=False,
        terminate_on_apogee=True,
    )

    # With Airbrakes
    sim_airbrakes = AirbrakesSimulation()
    sim_airbrakes.add_airbrakes(num_airbrakes=4, airbrakes_area=0.02)
    test_flight_airbrakes = Flight(
        rocket=sim_airbrakes.calisto,
        environment=sim_airbrakes.env,
        rail_length=5.2,
        inclination=85,
        heading=0,
        time_overshoot=False,
        terminate_on_apogee=True,
    )

    test_flight_airbrakes.rocket.evaluate_center_of_mass()
    test_flight_airbrakes.rocket.evaluate_center_of_pressure()
    print("Center of mass", test_flight_airbrakes.rocket.center_of_mass)
    print("Center of mass", test_flight_airbrakes.rocket.evaluate_center_of_pressure())

    print("Flight w/o Airbrakes Apogee:", test_flight_normal.apogee)
    print("Flight w Airbrakes Apogee:", test_flight_airbrakes.apogee)

    # time_list, deployment_level_list, drag_coefficient_list = [], [], []
    #
    # obs_vars = test_flight_airbrakes.get_controller_observed_variables()
    #
    # print(test_flight_airbrakes.apogee)
    #
    # for time, deployment_level, drag_coefficient in obs_vars:
    #     time_list.append(time)
    #     deployment_level_list.append(deployment_level)
    #     drag_coefficient_list.append(drag_coefficient)
    #
    # # Plot deployment level by time
    # plt.plot(time_list, deployment_level_list)
    # plt.xlabel("Time (s)")
    # plt.ylabel("Deployment Level")
    # plt.title("Deployment Level by Time")
    # plt.grid()
    # plt.show()
    #
    # # Plot drag coefficient by time
    # plt.plot(time_list, drag_coefficient_list)
    # plt.xlabel("Time (s)")
    # plt.ylabel("Drag Coefficient")
    # plt.title("Drag Coefficient by Time")
    # plt.grid()
    # plt.show()
