def create_contact_sphere_force_table(model, solution):
    """
    Create two tables of contact sphere forces and CoP location from the solution
    """
    model.initSystem()
    external_force_table = osim.TimeSeriesTableVec3()
    cop_table = osim.TimeSeriesTableVec3()
    
    states_trajectory = solution.exportToStatesTrajectory(model)
    num_states = states_trajectory.getSize()

    force_names = ['forceset/contactHeel',
                      'forceset/contactLateralRearfoot',
                      'forceset/contactLateralMidfoot',
                      'forceset/contactLateralToe',
                      'forceset/contactMedialToe',
                      'forceset/contactMedialMidfoot'
                      ]
    
    force_labels = ['heel',
                       'lat_rear',
                       'lat_mid', 
                       'lat_toe', 
                       'med_toe', 
                       'med_mid',
                       ]
    
    sphere_names = ['contactgeometryset/heel',
                       'contactgeometryset/lateralRearfoot',
                       'contactgeometryset/lateralMidfoot',
                       'contactgeometryset/lateralToe',
                       'contactgeometryset/medialToe',
                       'contactgeometryset/medialMidfoot',
                       ]
    
    previous_time = None

    for istate in range(num_states):
        state = states_trajectory.get(istate)
        # current_time = float(time_vector[istate])
        current_time = state.getTime()

        if previous_time is not None and abs(current_time - previous_time) < 1e-10:
            # Skip this state if the time is the same as the previous state
            continue

        # create state from the tabel data
        state = model.initSystem()


        # state = states_table.get(istate)
        model.realizeVelocity(state)
        # model.realizeDynamics(state)
        
        row = osim.RowVectorVec3(3*2*len(force_names))
        cop_row = osim.RowVectorVec3((2))
        labels = osim.StdVectorString()
        for iside, side in enumerate(['_r', '_l']):

            cop = osim.Vec3(0)
            cop_torque_sum_x = 0
            cop_torque_sum_z = 0
            cop_force_sum = 0

            offset = 3 * iside * len(sphere_names)
            zipped = zip(force_names, force_labels, sphere_names)
            
            for i, (force_name, force_label, sphere_name) in enumerate(zipped):

                force = osim.Vec3(0)
                torque = osim.Vec3(0)

                force_obj = osim.Force.safeDownCast(model.getComponent(f'{force_name}{side}'))
                force_vals = force_obj.getRecordValues(state)

                force[0] = force_vals.get(0)
                force[1] = force_vals.get(1)
                force[2] = force_vals.get(2)
                torque[0] = force_vals.get(3)
                torque[1] = force_vals.get(4)
                torque[2] = force_vals.get(5)

                sphere = osim.ContactSphere.safeDownCast(model.getComponent(f'{sphere_name}{side}'))
                frame = sphere.getFrame()
                position = frame.getPositionInGround(state)
                location = sphere.get_location()
                location_in_ground = frame.expressVectorInGround(state, location)

                position[0] = position[0] + location_in_ground[0]
                position[1] = position[1] + location_in_ground[1]
                position[2] = position[2] + location_in_ground[2]
                
                row[3 * i + offset] = force
                row[3 * i + 1 + offset] = position
                row[3 * i + 2 + offset] = torque

                cop_force_sum += force[1]
                cop_torque_sum_x += force[1] * position[0]
                cop_torque_sum_z += force[1] * position[2]
                
                for suffix in ['_force_v', '_force_p', '_torque_']:
                    labels.append(f'{force_label}{side}{suffix}')

            if np.abs(cop_force_sum) > 0.001:
                cop[0] = cop_torque_sum_x / cop_force_sum
                cop[2] = cop_torque_sum_z / cop_force_sum

            cop_row[iside] = cop

            external_force_table.appendRow(state.getTime(), row)
            cop_table.appendRow(state.getTime(), cop_row)

        external_force_table.setcolumnLabels(labels)

        labels = osim.StdVectorString()
        labels.append('2_ground_force_p')
        labels.append('1_ground_force_p')
        cop_table.setColumnLabels(labels)

        suffixes = osim.StdVectorString()
        suffixes.append('x')
        suffixes.append('y')
        suffixes.append('z')

        return external_force_table.flatten(suffixes), cop_table

