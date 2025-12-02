# -*- coding: utf-8 -*-
"""

@author: 
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Build of functions to assist with running simulations
    
"""

import pandas as pd
import numpy as np
import os
import opensim as osim

# %% Function to add set of torque actuators to model

def addTorqueActuators(osimModel = None,
                       optForces = None,
					   controlLimits = None, 
                       reserve_values = None):
    
    """
    
    Convenience function for adding series of torque actuators to model
    
    Input:    osimModel - OpenSim model object for use
              optForces - dict of coordinates and their associated optimal forces to add
			  controlLimits - dict of coordinates and their associated max/min control limits to add
              reserve_values - single value to override prescribed dictionary values in optForces arg, sets this value as opt force for all coordinate actuators/reserves (optional)

    Output:   osimModel - updated torque driven model
                  
    """
    
    #Check inputs
    if osimModel is None or optForces is None:
            raise ValueError('All inputs for this function are required!')

    #Intialise model system
    osimModel.initSystem()

    if reserve_values is not None:
        # generate new dict where all values are single reserve
        single_optForces = {coord: reserve_values for coord in optForces}
        optForces = single_optForces
    
    #Get coordinate list
    coordinatesList = list(optForces.keys())
    
    #Get coordinate set
    coordSet = osimModel.getCoordinateSet()
    
    #Loop through coordinates and add actuators
    for coordinate in coordinatesList:
        #Create actuator
        actu = osim.CoordinateActuator()
        #Set name
        actu.setName(f'{coordinate}_actuator')
        #Set coordinate
        actu.setCoordinate(coordSet.get(coordinate))
        #Set optimal force
        actu.setOptimalForce(optForces[coordinate])
        #Set min and max control
        actu.setMinControl(controlLimits[coordinate]*-1)
        actu.setMaxControl(controlLimits[coordinate]*1)
        #Append to model force set
        osimModel.updForceSet().cloneAndAppend(actu)
    
    #Finalise model connections
    osimModel.finalizeConnections()
    
    #Return model
    return osimModel

# %% Function to convert IK coordinates to states

def kinematicsToStates(kinematicsFileName=None, osimModelFileName=None,
                       outputFileName=None,
                       inDegrees=True, outDegrees=False,
                       filtFreq=None,
                       excludeColumns=None,
                       include_speeds=False,
                       include_accelerations=False):
    """
    Enhanced convenience function for converting IK results to a states storage using TimeSeriesTable.
    Can optionally compute and include angular velocities and accelerations.
    
    Parameters:
    -----------
    kinematicsFileName : str
        File containing kinematic data. Header should only be coordinates name, rather than path to state
    osimModelFileName : str  
        OpenSim model filename that corresponds to kinematic data
    outputFileName : str, optional
        Filename to output to (defaults to 'coordinates.sto')
    inDegrees : bool, optional
        Set to True if kinematics file is in degrees (defaults to True)
    outDegrees : bool, optional
        Set to True if desired output is in degrees (defaults to False)
    filtFreq : float, optional
        Lowpass filter frequency for kinematics (optional)
    excludeColumns : list, optional
        List of column name patterns to exclude from degree/radian conversion 
        (defaults to ['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz'])
    include_speeds : bool, optional
        Set to True to compute and include angular velocities (defaults to False)
    include_accelerations : bool, optional
        Set to True to compute and include angular accelerations (defaults to False)
    """
    
    if kinematicsFileName is None:
        raise ValueError('Filename for kinematics is required')
    if osimModelFileName is None:
        raise ValueError('OpenSim model filename is required')
    
    # Load the model to get coordinate information
    kinematicModel = osim.Model(osimModelFileName)
    kinematicModel.initSystem()
    coords_set = kinematicModel.getCoordinateSet()
    
    
    # Set default exclude columns if not provided
    if excludeColumns is None:
        excludeColumns = ['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz']

    print(f"Enhanced kinematicsToStates processing: {kinematicsFileName}")
    if include_speeds:
        print("  - Will compute angular velocities")
    if include_accelerations:
        print("  - Will compute angular accelerations")

    # Load the kinematic data using TimeSeriesTable
    kinematicsTable = osim.TimeSeriesTable(kinematicsFileName)

    # Apply filtering if requested
    if filtFreq is not None:
        print(f"  - Applying lowpass filter at {filtFreq} Hz")
        # Create a TableProcessor for filtering
        tableProcessor = osim.TableProcessor(kinematicsFileName)
        tableProcessor.append(osim.TabOpLowPassFilter(filtFreq))
        kinematicsTable = tableProcessor.processAndConvertToRadians(kinematicModel)

    # Get time data
    time_data = kinematicsTable.getIndependentColumn()
    time_array = np.array([time_data[i] for i in range(len(time_data))])
    
    # Create output table
    time_vector = osim.StdVectorDouble()
    for t in time_array:
        time_vector.push_back(t)
    outputTable = osim.TimeSeriesTable(time_vector)
    
    # Get column labels and process coordinate information
    columnLabels = kinematicsTable.getColumnLabels()
    coordinate_info = []  # Store (coord_name, coord_object, full_path, data) for derivative computation
    
    # Process each column and rename to full paths
    for i in range(len(columnLabels)):
        currAngle = columnLabels[i]
        print(f"processing {len(columnLabels)} columns from input file..")
        
        # Handle both moco solution full paths and standard ik names of just coordinate name
        if '/jointset/' in currAngle and '/value' in currAngle:
            # for full path moco inputs
            path_parts = currAngle.split('/')
            if len(path_parts) >= 4:
                coord_name = path_parts[-2]
                fullPath = currAngle
            else:
                print(f" Can't parse coordinate path: {currAngle}")
                continue
        else:
            # for standard ik column labels
            coord_name = currAngle
            fullPath = None # (gets construvted after finding coordinate)

        # Get the original data for given column
        try:
            columnData = kinematicsTable.getDependentColumn(currAngle)
            data = columnData.to_numpy()
        except Exception as e:
            print(f" Can't get data for column '{currAngle}': {e}")
            continue

        # Try getting the full path to coordinate
        try:
            # Look for coordinate in model
            coordinate = coords_set.get(coord_name)

            # coordinate = kinematicModel.updCoordinateSet().get(currAngle)
            if fullPath is None:
                fullPath = coordinate.getAbsolutePathString() + '/value'
                print(f" Found coord: {coord_name} -> {fullPath}")

            # Check if this column should be excluded from conversion
            shouldExclude = any(excl in currAngle for excl in excludeColumns)

            # Get column data
            # columnData = kinematicsTable.getDependentColumn(currAngle).to_numpy()
            convertedData = data.copy()

            if not shouldExclude:
                # Handle unit conversions for position data
                if inDegrees and not outDegrees:
                    # Convert degrees to radians
                    convertedData = data * (np.pi / 180.0)
                    print(f"  - Converted {currAngle} from degrees to radians")
                
                elif not inDegrees and outDegrees:
                    # Convert radians to degrees
                    convertedData = data * (180.0 / np.pi)
                    print(f"  - Converted {currAngle} from radians to degrees")
                else:    
                    convertedData = data
            else:
                print(f"  - Excluding {currAngle} from radian to degree conversion")
                
            # Add to output table
            converted_dataVector = osim.Vector(convertedData.tolist())
            outputTable.appendColumn(fullPath, converted_dataVector)
            
            # Store coordinate info for derivative computation
            coordinate_info.append((currAngle, coordinate, fullPath, convertedData))
            print(f"  - Found coordinate: {currAngle} -> {fullPath}")
            
        except Exception as e:
            # Print out that current column isn't a coordinate
            print(f"  - Can't fincd coord: {currAngle} in model: {e}")
            print(f" (Original Column: {currAngle})")
            continue
    print(f" successfully processed {len(coordinate_info)} coordinates")
    
    # Compute derivatives if requested
    if include_speeds or include_accelerations:
        print(f"  - Computing derivatives for {len(coordinate_info)} coordinates...")
        
        # Calculate time step (assuming uniform sampling)
        dt = time_array[1] - time_array[0] if len(time_array) > 1 else 1.0
        
        # Compute speed derivatives
        if include_speeds:
            print(f"  - Adding {len(coordinate_info)} speed columns to output")
            for coord_name, coord_obj, full_path, data in coordinate_info:
                # Calculate velocity using numpy gradient
                velocity_data = np.gradient(data, dt)
                
                # Apply unit conversion to velocities if needed
                if inDegrees and not outDegrees:
                    shouldExclude = any(excl in coord_name for excl in excludeColumns)
                    if not shouldExclude:
                        velocity_data = velocity_data * (np.pi / 180.0)
                elif not inDegrees and outDegrees:
                    shouldExclude = any(excl in coord_name for excl in excludeColumns)
                    if not shouldExclude:
                        velocity_data = velocity_data * (180.0 / np.pi)
                
                # Create speed column name and add to table
                speed_path = full_path.replace('/value', '/speed')
                velocityVector = osim.Vector(velocity_data.tolist())
                outputTable.appendColumn(speed_path, velocityVector)
        
        # Compute acceleration derivatives
        if include_accelerations:
            print(f"  - Adding {len(coordinate_info)} acceleration columns to output")
            for coord_name, coord_obj, full_path, data in coordinate_info:
                # Calculate velocity first, then acceleration
                velocity_data = np.gradient(data, dt)
                acceleration_data = np.gradient(velocity_data, dt)
                
                # Apply unit conversion to accelerations if needed
                if inDegrees and not outDegrees:
                    shouldExclude = any(excl in coord_name for excl in excludeColumns)
                    if not shouldExclude:
                        acceleration_data = acceleration_data * (np.pi / 180.0)
                elif not inDegrees and outDegrees:
                    shouldExclude = any(excl in coord_name for excl in excludeColumns)
                    if not shouldExclude:
                        acceleration_data = acceleration_data * (180.0 / np.pi)
                
                # Create acceleration column name and add to table
                accel_path = full_path.replace('/value', '/acceleration')
                accelVector = osim.Vector(acceleration_data.tolist())
                outputTable.appendColumn(accel_path, accelVector)
    
    # Set table metadata
    if outDegrees:
        outputTable.addTableMetaDataString('inDegrees', 'yes')
    else:
        outputTable.addTableMetaDataString('inDegrees', 'no')
    
    # Write the output table to file
    osim.STOFileAdapter.write(outputTable, outputFileName)
    
    print(f"States file written to: {outputFileName}")
    print(f"  - Total coordinates processed: {len(coordinate_info)}")
    if include_speeds:
        print(f"  - Speed columns added: {len(coordinate_info)}")
    if include_accelerations:
        print(f"  - Acceleration columns added: {len(coordinate_info)}")


    
# %% Function to convert states kinematics back to standard kinematics
 
def kinematicsToStatesForMoco(mocoSolutionFileName, outputFileName='kinematics.sto', 
                              inDegrees=False, outDegrees=True, verbose=True):
    """
    Reads a Moco solution .sto file and converts coordinate values from radians to degrees.
    Only processes columns that:
    - Start with '/jointset/'
    - End with '/value'
    - Are NOT 'pelvis_tx/value', 'pelvis_ty/value', or 'pelvis_tz/value'
    
    Parameters:
    -----------
    mocoSolutionFileName : str
        Path to the Moco solution .sto file (typically in radians)
    outputFileName : str
        Path to save the converted kinematics file (in degrees)
    inDegrees : bool
        Set to True if input is already in degrees (default: False)
    outDegrees : bool
        Set to True to convert output to degrees (default: True)
    verbose : bool
        Print processing information (default: True)
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing the converted kinematics data
    """
    
    if verbose:
        print(f"Reading Moco solution from: {mocoSolutionFileName}")
    
    # Read the .sto file line by line
    with open(mocoSolutionFileName, 'r') as f:
        lines = f.readlines()
    
    # Find the 'endheader' line
    header_end_idx = None
    header_lines = []
    
    for i, line in enumerate(lines):
        if 'endheader' in line.lower():
            header_end_idx = i
            break
        header_lines.append(line.rstrip())
    
    if header_end_idx is None:
        raise ValueError(f"Could not find 'endheader' in file: {mocoSolutionFileName}")
    
    if verbose:
        print(f"Found 'endheader' at line {header_end_idx}")
    
    # The line after 'endheader' contains column headers
    column_headers_line = lines[header_end_idx + 1].strip()
    column_headers = column_headers_line.split('\t')
    
    if verbose:
        print(f"Found {len(column_headers)} columns")
    
    # Parse the data starting from two lines after 'endheader'
    data_start_idx = header_end_idx + 2
    data_lines = lines[data_start_idx:]
    
    # Parse data into a list of lists
    data_rows = []
    for line in data_lines:
        line = line.strip()
        if line:  # Skip empty lines
            try:
                values = [float(x) for x in line.split('\t')]
                data_rows.append(values)
            except ValueError:
                if verbose:
                    print(f"Warning: Could not parse line: {line[:50]}...")
                continue
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, columns=column_headers)
    
    if verbose:
        print(f"Loaded {len(df)} rows of data")
    
    # Define columns to exclude from conversion (pelvis translations)
    exclude_coords = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz']
    
    # Identify columns to convert
    columns_to_convert = []
    
    for col in column_headers:
        # Check if column matches criteria:
        # 1. Starts with '/jointset/'
        # 2. Ends with '/value'
        # 3. Not a pelvis translation
        
        if col.startswith('/jointset/') and col.endswith('/value'):
            # Check if it's not a pelvis translation
            should_exclude = any(excl in col for excl in exclude_coords)
            
            if not should_exclude:
                columns_to_convert.append(col)
    
    if verbose:
        print(f"\nColumns to convert from radians to degrees: {len(columns_to_convert)}")
        for col in columns_to_convert:
            print(f"  - {col}")
    
    # Perform conversion if needed
    if not inDegrees and outDegrees:
        if verbose:
            print("\nConverting from radians to degrees...")
        
        for col in columns_to_convert:
            # Convert radians to degrees
            df[col] = df[col] * (180.0 / np.pi)
            if verbose:
                print(f"  Converted: {col}")
    
    elif inDegrees and not outDegrees:
        if verbose:
            print("\nConverting from degrees to radians...")
        
        for col in columns_to_convert:
            # Convert degrees to radians
            df[col] = df[col] * (np.pi / 180.0)
            if verbose:
                print(f"  Converted: {col}")
    
    else:
        if verbose:
            print("\nNo unit conversion needed (inDegrees={}, outDegrees={})".format(
                inDegrees, outDegrees))
    
    # Write output file
    if verbose:
        print(f"\nWriting converted file to: {outputFileName}")
    
    with open(outputFileName, 'w') as f:
        # Write header
        for line in header_lines:
            # Update inDegrees metadata if present
            if 'inDegrees' in line.lower():
                if outDegrees:
                    f.write('inDegrees=yes\n')
                else:
                    f.write('inDegrees=no\n')
            else:
                f.write(line + '\n')
        
        # Add inDegrees metadata if not present in original
        if not any('inDegrees' in line.lower() for line in header_lines):
            if outDegrees:
                f.write('inDegrees=yes\n')
            else:
                f.write('inDegrees=no\n')
        
        # Write endheader
        f.write('endheader\n')
        
        # Write column headers
        f.write('\t'.join(column_headers) + '\n')
        
        # Write data
        df.to_csv(f, sep='\t', index=False, header=False, float_format='%.8f')
    
    if verbose:
        print(f"\n✅ Conversion completed successfully")
        print(f"Output file: {outputFileName}")
        print(f"Rows: {len(df)}, Columns: {len(df.columns)}")
    
    return df
 
def compute_sto_stdevs(stoFileName=None, outputCsvFileName='standard_deviations.csv'):
    """
    Calculate standard deviations for each column in a .sto file and write to CSV.
    
    This function reads a .sto file (typically output from kinematicsToStates),
    calculates the standard deviation for each column (excluding time), and
    outputs the results to a CSV file with column labels and their corresponding
    standard deviation values.
    
    Input:    stoFileName - path to the .sto file to analyze
              outputCsvFileName - optional filename for CSV output (defaults to 'standard_deviations.csv')
              
    Output:   CSV file with two columns: 'Column_Label' and 'Standard_Deviation'
    """
    
    if stoFileName is None:
        raise ValueError('Filename for .sto file is required')
    
    # Check if the .sto file exists
    if not os.path.exists(stoFileName):
        raise FileNotFoundError(f'File {stoFileName} not found')
    
    try:
        # Load the .sto file using OpenSim
        # Read the .sto file as a TimeSeriesTable
        stoTable = osim.TimeSeriesTable(stoFileName)
        
        # Get column labels
        columnLabels = stoTable.getColumnLabels()
        
        # Initialize lists to store results
        labels = []
        std_values = []
        
        # Handle different OpenSim versions - columnLabels might be tuple or OpenSim object
        if hasattr(columnLabels, 'getSize'):
            # OpenSim object with getSize() method
            num_columns = columnLabels.getSize()
            get_column_name = lambda i: columnLabels.get(i)
            # for i in range(num_columns):
            #     current_label = get_column_name(i)
            #     if '/value' in current_label:
            #         print(f"Processing column: {current_label}")
            #         pass
        else:
            # Python tuple/list
            num_columns = len(columnLabels)
            get_column_name = lambda i: columnLabels[i]
            # for i in range(num_columns):
            #     current_label = get_column_name(i)
            #     if '/value' in current_label:
            #         print(f"Processing column: {current_label}")
            #         pass
        
        # Calculate standard deviation for each column
        for i in range(num_columns):
            columnName = get_column_name(i)
            
            # Extract text between third and fourth '/' characters
            # Split by '/' and get the appropriate part
            if '/value' not in columnName: continue
            parts = str(columnName).split('/')
            
            if len(parts) >= 4:
                # Get the text between 3rd and 4th '/' (index 3, since split creates parts)
                extracted_name = parts[3]
            else:
                # If there aren't enough '/' characters, use the original name
                extracted_name = str(columnName)
            
            # Get the column data
            columnData = stoTable.getDependentColumn(columnName)
            
            # Convert to numpy array for standard deviation calculation
            # Handle different OpenSim data types
            if hasattr(columnData, 'size') and hasattr(columnData, 'get'):
                # OpenSim Vector object with get() method
                dataArray = np.array([columnData.get(j) for j in range(columnData.size())])
            elif hasattr(columnData, '__len__') and hasattr(columnData, '__getitem__'):
                # VectorView or similar - can be accessed like a Python sequence
                dataArray = np.array([columnData[j] for j in range(len(columnData))])
            elif hasattr(columnData, 'to_numpy'):
                # Some OpenSim objects have to_numpy() method
                dataArray = columnData.to_numpy()
            else:
                # Try direct conversion to numpy array
                try:
                    dataArray = np.array(columnData)
                except:
                    # Last resort - try to convert to list first
                    dataArray = np.array(list(columnData))
            
            # Calculate standard deviation
            stdDev = np.std(dataArray, ddof=1)  # Using sample standard deviation (ddof=1)
            
            # Store results using the extracted name
            labels.append(extracted_name)
            std_values.append(stdDev)
            
            print(f'Standard deviation for {extracted_name} (from {columnName}): {stdDev:.6f}')
        
    except ImportError:
        # Fallback method if OpenSim is not available - parse manually
        print("OpenSim not available, using manual parsing method...")
        
        # Read file line by line to find 'endheader'
        with open(stoFileName, 'r') as f:
            lines = f.readlines()
        
        # Find the line with 'endheader'
        endheader_idx = None
        for i, line in enumerate(lines):
            if 'endheader' in line:
                endheader_idx = i
                break
        
        if endheader_idx is None:
            raise ValueError(f"No 'endheader' found in file: {stoFileName}")
        
        # Headers are in the line after 'endheader'
        header_line = lines[endheader_idx + 1]
        
        # Data starts two lines after 'endheader'
        data_start_idx = endheader_idx + 2
        
        # Parse headers (tab-delimited)
        headers = [h.strip() for h in header_line.split('\t') if h.strip()]
        
        # Create a dictionary to store the data
        data_dict = {header: [] for header in headers}
        
        # Parse the data
        for i in range(data_start_idx, len(lines)):
            line = lines[i].strip()
            if not line:  # Skip empty lines
                continue
            
            values = line.split('\t')
            
            # Make sure we have the right number of values
            if len(values) >= len(headers):
                for j, header in enumerate(headers):
                    try:
                        data_dict[header].append(float(values[j]))
                    except (ValueError, IndexError):
                        # Handle non-numeric data or missing values
                        pass
        
        # Convert to pandas DataFrame
        df = pd.DataFrame(data_dict)
        
        # Initialize lists to store results
        labels = []
        std_values = []
        
        # Calculate standard deviation for each column (excluding time if present)
        for column in df.columns:
            if column.lower() != 'time':  # Skip time column
                # Extract text between third and fourth '/' characters
                parts = str(column).split('/')
                if len(parts) >= 4:
                    # Get the text between 3rd and 4th '/' (index 3, since split creates parts)
                    extracted_name = parts[3]
                else:
                    # If there aren't enough '/' characters, use the original name
                    extracted_name = str(column)
                
                # Calculate standard deviation
                stdDev = df[column].std(ddof=1)  # Using sample standard deviation
                
                # Store results using the extracted name
                labels.append(extracted_name)
                std_values.append(stdDev)
                
                print(f'Standard deviation for {extracted_name} (from {column}): {stdDev:.6f}')
    
    # Clean up labels and values - remove any invalid entries
    clean_labels = []
    clean_std_values = []
    
    for label, std_val in zip(labels, std_values):
        # Skip if label is None, empty, or std_val is NaN
        if label is not None and str(label).strip() != '' and not np.isnan(std_val):
            clean_labels.append(str(label).strip())
            clean_std_values.append(std_val)
        else:
            print(f"Skipping invalid entry: label='{label}', std_dev={std_val}")
    
    # Create results DataFrame
    results_df = pd.DataFrame({
        '': clean_labels,
        '0': clean_std_values
    })
    
    # Ensure no duplicate labels
    # results_df = results_df.drop_duplicates(subset=['Column_Label'])
    
    # Write to CSV file
    results_df.to_csv(outputCsvFileName, index=False)
    
    print(f'\nStandard deviations written to: {outputCsvFileName}')
    print(f'Total columns analyzed: {len(labels)}')
    
    return results_df

def generate_guess_file(kinematicsFileName=None, osimModelFileName=None,
                        outputFileName='guess_file.sto',
                        muscle_activation_bounds=[0.01, 1.0],
                        normalized_tendon_force_bounds=[0.01, 1.8],
                        torque_actuator_bounds=[-1.0, 1.0],
                        inDegrees=True, filtFreq=None):
    """
    Generate a guess file from kinematic data with specific requirements:
    - Time column
    - Coordinate values (converted to radians) with full paths
    - First derivative (speed) values for each coordinate
    - Muscle activation columns for ALL muscles (default values)
    - Normalized tendon force columns ONLY for specific ankle plantarflexors
    - Torque actuator activation and excitation columns ONLY for upper body coordinates
    
    Input:    kinematicsFileName - file containing kinematic data (.sto file)
              osimModelFileName - opensim model filename that corresponds to kinematic data
              outputFileName - filename to output to (defaults to 'guess_file.sto')
              muscle_activation_bounds - list with [lower_bound, upper_bound] for muscle activations
              normalized_tendon_force_bounds - list with [lower_bound, upper_bound] for specific tendon forces
              torque_actuator_bounds - list with [lower_bound, upper_bound] for upper body torque actuators
              inDegrees - set to true if kinematics file is in degrees (defaults to True)
              filtFreq - lowpass filter frequency for kinematics
    
    Output:   Writes .sto file with selective guess data
    """
    
    if kinematicsFileName is None:
        raise ValueError('Filename for kinematics is required')
    if osimModelFileName is None:
        raise ValueError('OpenSim model filename is required')
    
    # Load the model to get coordinate and muscle information
    model = osim.Model(osimModelFileName)
    model.initSystem()
    
    # Use Storage class for reliable data handling
    kinematicsStorage = osim.Storage(kinematicsFileName)
    
    # Apply filtering if requested
    if filtFreq is not None:
        kinematicsStorage.lowpassFIR(4, filtFreq)
    
    # Get coordinate set for coordinate validation
    coordSet = model.getCoordinateSet()
    
    # Get time data
    timeArray = osim.ArrayDouble()
    kinematicsStorage.getTimeColumn(timeArray)
    times = [timeArray.get(i) for i in range(timeArray.getSize())]
    numRows = len(times)
    
    # Create output table with time vector
    timeVector = osim.StdVectorDouble()
    for time_val in times:
        timeVector.push_back(time_val)
    outputTable = osim.TimeSeriesTable(timeVector)
    
    # Get original column labels
    originalLabels = kinematicsStorage.getColumnLabels()
    coordinateColumns = []  # Track coordinate info for derivative calculation
    excludeColumns = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz']
    
    print(f"Processing {originalLabels.getSize()} columns from input file...")
    
    # Process each original column (coordinates)
    for i in range(originalLabels.getSize()):
        if originalLabels.get(i) == 'time':
            continue
            
        origName = originalLabels.get(i)
        
        # Get the data for this column
        dataArray = osim.ArrayDouble()
        kinematicsStorage.getDataColumn(origName, dataArray)
        data = [dataArray.get(j) for j in range(dataArray.getSize())]
        
        # Try to find this coordinate in the model
        try:
            coord = coordSet.get(origName)
            fullPath = coord.getAbsolutePathString() + '/value'
            
            print(f"Found coordinate: {origName} -> {fullPath}")
            
            # Convert degrees to radians if needed (excluding pelvis translations)
            if inDegrees and origName not in excludeColumns:
                # Convert each value from degrees to radians
                convertedData = [val * (np.pi / 180.0) for val in data]
                # Convert to OpenSim Vector
                convertedVector = osim.Vector(convertedData)
                outputTable.appendColumn(fullPath, convertedVector)
                coordinateColumns.append((origName, fullPath, convertedData))
                print(f"  Converted {origName} from degrees to radians")
            else:
                # Keep original values (for pelvis translations or if already in radians)
                originalVector = osim.Vector(data)
                outputTable.appendColumn(fullPath, originalVector)
                coordinateColumns.append((origName, fullPath, data))
                print(f"  Kept {origName} in original units")
                
        except:
            # If not found as direct coordinate name, try alternative approaches
            coordinate_found = False
            
            # Try to find coordinate by iterating through all coordinates in the model
            for j in range(coordSet.getSize()):
                coord = coordSet.get(j)
                coordName = coord.getName()
                
                # Check if this might be a match (case-insensitive, handle underscores)
                if (coordName.lower() == origName.lower() or 
                    coordName.lower().replace('_', '') == origName.lower().replace('_', '') or
                    origName.lower() in coordName.lower() or
                    coordName.lower() in origName.lower()):
                    
                    fullPath = coord.getAbsolutePathString() + '/value'
                    print(f"Found coordinate (alternative): {origName} -> {coordName} -> {fullPath}")
                    
                    if inDegrees and origName not in excludeColumns:
                        convertedData = [val * (np.pi / 180.0) for val in data]
                        convertedVector = osim.Vector(convertedData)
                        outputTable.appendColumn(fullPath, convertedVector)
                        coordinateColumns.append((origName, fullPath, convertedData))
                        print(f"  Converted {origName} from degrees to radians")
                    else:
                        originalVector = osim.Vector(data)
                        outputTable.appendColumn(fullPath, originalVector)
                        coordinateColumns.append((origName, fullPath, data))
                        print(f"  Kept {origName} in original units")
                    
                    coordinate_found = True
                    break
            
            if not coordinate_found:
                print(f"Warning: {origName} not found as a coordinate in the model...skipping...")
    
    # Calculate and add speed derivatives for coordinates
    for origName, fullPath, data in coordinateColumns:
        # Calculate derivatives using central differences
        derivatives = []
        for i in range(len(data)):
            if i == 0:
                # Forward difference for first point
                deriv = (data[i+1] - data[i]) / (times[i+1] - times[i])
            elif i == len(data) - 1:
                # Backward difference for last point
                deriv = (data[i] - data[i-1]) / (times[i] - times[i-1])
            else:
                # Central difference for interior points
                deriv = (data[i+1] - data[i-1]) / (times[i+1] - times[i-1])
            derivatives.append(deriv)
        
        # Create speed column name and add to table
        speedPath = fullPath.replace('/value', '/speed')
        derivativesVector = osim.Vector(derivatives)
        outputTable.appendColumn(speedPath, derivativesVector)
        print(f"Added speed derivative: {speedPath}")
    
    # Add muscle activation columns for ALL muscles
    muscleSet = model.getMuscles()
    activationValue = (muscle_activation_bounds[0] + muscle_activation_bounds[1]) / 2.0
    
    print(f"Adding muscle activations for {muscleSet.getSize()} muscles...")
    for i in range(muscleSet.getSize()):
        muscle = muscleSet.get(i)
        muscleName = muscle.getName()
        activationPath = f'/forceset/{muscleName}/activation'
        
        # Create constant activation values
        activationData = [activationValue] * numRows
        activationVector = osim.Vector(activationData)
        outputTable.appendColumn(activationPath, activationVector)
    
    # Add normalized tendon force columns ONLY for specific ankle plantarflexors
    tendonForceValue = (normalized_tendon_force_bounds[0] + normalized_tendon_force_bounds[1]) / 2.0
    plantarflexor_muscles = ['gas_med_l', 'gas_lat_l', 'soleus_l', 'gas_med_r', 'gas_lat_r', 'soleus_r']
    
    print(f"Adding normalized tendon forces for specific plantarflexors...")
    tendon_force_added = 0
    for i in range(muscleSet.getSize()):
        muscle = muscleSet.get(i)
        muscleName = muscle.getName()
        
        # Only add tendon force for specific plantarflexor muscles
        if muscleName in plantarflexor_muscles:
            tendonForcePath = f'/forceset/{muscleName}/normalized_tendon_force'
            
            # Create constant tendon force values
            tendonForceData = [tendonForceValue] * numRows
            tendonForceVector = osim.Vector(tendonForceData)
            outputTable.appendColumn(tendonForcePath, tendonForceVector)
            print(f"  Added tendon force for: {muscleName}")
            tendon_force_added += 1
    
    print(f"Total tendon force columns added: {tendon_force_added}")
    
    # Add torque actuator columns ONLY for upper body coordinates
    upper_body_coordinates = [
        'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
        'arm_flex_r', 'arm_add_r', 'arm_rot_r', 'elbow_flex_r', 'pro_sup_r',
        'arm_flex_l', 'arm_add_l', 'arm_rot_l', 'elbow_flex_l', 'pro_sup_l'
    ]
    
    forceSet = model.getForceSet()
    torqueActivationValue = (torque_actuator_bounds[0] + torque_actuator_bounds[1]) / 2.0
    torqueExcitationValue = (torque_actuator_bounds[0] + torque_actuator_bounds[1]) / 2.0
    
    print(f"Adding torque actuator columns for upper body coordinates...")
    torque_actuators_added = 0
    
    # Look for torque actuators matching upper body coordinates
    for i in range(forceSet.getSize()):
        force = forceSet.get(i)
        forceName = force.getName()
        
        # Check if it's a torque actuator for upper body coordinates
        if 'torque_' in forceName.lower():
            # Extract coordinate name from torque actuator name (e.g., 'torque_arm_flex_r' -> 'arm_flex_r')
            coord_name = forceName.lower().replace('torque_', '')
            
            if coord_name in upper_body_coordinates:
                # Add activation column
                activationPath = f'/forceset/{forceName}/activation'
                activationData = [torqueActivationValue] * numRows
                activationVector = osim.Vector(activationData)
                outputTable.appendColumn(activationPath, activationVector)
                
                # Add excitation/control column
                excitationPath = f'/forceset/{forceName}'
                excitationData = [torqueExcitationValue] * numRows
                excitationVector = osim.Vector(excitationData)
                outputTable.appendColumn(excitationPath, excitationVector)
                
                print(f"  Added torque actuator: {forceName}")
                torque_actuators_added += 1
    
    print(f"Total torque actuator pairs added: {torque_actuators_added}")
    
    # Set table metadata
    outputTable.addTableMetaDataString('inDegrees', 'no')  # All angles now in radians
        # Add required metadata for MocoTrajectory
    total_columns = len(coordinateColumns) * 2 + muscleSet.getSize() + tendon_force_added + torque_actuators_added
    outputTable.addTableMetaDataString('num_states', str(total_columns))
    outputTable.addTableMetaDataString('num_controls', '0')  # No controls in states file
    outputTable.addTableMetaDataString('num_derivatives', '0')  # No derivatives tracking
    outputTable.addTableMetaDataString('num_parameters', '0')  # No parameters
    outputTable.addTableMetaDataString('num_multipliers', '0')  # No multipliers
    outputTable.addTableMetaDataString('num_slacks', '0')  # No multipliers
    outputTable.addTableMetaDataString('num_input_controls', '0')  # No input controls

    # Optional: Add solver metadata for completeness
    outputTable.addTableMetaDataString('solver', 'generated_guess')
    outputTable.addTableMetaDataString('write_version', '4.0')
    outputTable.addTableMetaDataString('problem_type', 'tracking')
    # Write the table to file
    osim.STOFileAdapter.write(outputTable, outputFileName)
    
    print(f"\n✅ Generated selective guess file: {outputFileName}")
    print(f"=== SUMMARY ===")
    print(f"  - Input coordinates processed: {len(coordinateColumns)}")
    print(f"  - Speed derivatives added: {len(coordinateColumns)}")
    print(f"  - Muscle activations (ALL muscles): {muscleSet.getSize()} (default: {activationValue:.3f})")
    print(f"  - Normalized tendon forces (plantarflexors only): {tendon_force_added} (default: {tendonForceValue:.3f})")
    print(f"  - Upper body torque actuators: {torque_actuators_added} pairs (activation: {torqueActivationValue:.3f}, excitation: {torqueExcitationValue:.3f})")
    print(f"  - Total time points: {numRows}")
    print(f"  - Angle units: {'Converted to radians' if inDegrees else 'Already in radians'}")
    print(f"\n🚫 EXCLUDED:")
    print(f"  - No reserve actuators")
    print(f"  - No muscle excitations/controls")
    print(f"  - Tendon forces only for specific plantarflexors")
    print(f"  - Torque actuators only for upper body coordinates")
    
def statesToKinematics(statesFileName = None, osimModelFileName = None,
                       outputFileName = 'coordinates.sto',
                       inDegrees = False, outDegrees = True):
    
    # Convenience function for converting IK results to a states storage.
    #
    # Input:    statesFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
    #           osimModelFileName - opensim model filename that corresponds to kinematic data
    #           outputFileName - optional filename to output to (defaults to coordinates.sto)
    #           inDegrees - set to true if kinematics file is in degrees (defaults to False)
    #           outDegrees - set to true if desired output is in degrees (defaults to True)

    if statesFileName is None:
        raise ValueError('Filename for states is required')

    excludeColumns = ['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz']
    #Load in the states data as a table
    statesTable = osim.TimeSeriesTable(statesFileName)
    
    # #Create a copy of the states data to alter the column labels in
    kinematicsStorage = osim.Storage(statesFileName)
    
    #Get the column headers for the states file
    stateNames = statesTable.getColumnLabels()
    
    kinematicModel = osim.Model(osimModelFileName)

    #Loop through the column names and identify which columns need to be removed
    #given they aren't kinematic values
    for state in stateNames:
        if '/value' not in state:
            statesTable.removeColumn(state)
    
    #Create new set of column labels that removes all but the coordinate
    #This is based off the fact that kinematic states are presented as:
    # /jointset/joint_name/coordinate_name/value
    
    #Get the selected state names and create a string array to store names in
    selectStateNames = statesTable.getColumnLabels()
    newStateNames = osim.StdVectorString()
    
    #Loop through state names to alter
    for state in selectStateNames:
        #Split the string by the / and get the 4th output
        #Append this to the vector string
        newStateNames.append(state.split('/')[2])
        
    #Set new names in table
    statesTable.setColumnLabels(newStateNames)

    #     #Appropriately set output in degrees or radians
    # if inDegrees and not outDegrees:
    #     #Convert degrees values to radians for selected columns only
    #     #Get all the column labels to identify which to convert
    #     columnLabels = statesStorage.getColumnLabels()
        
    #     print(f"Converting degrees to radians, excluding: {excludeColumns}")
        
    #     # Manually convert degrees to radians for all rows and columns
    #     for iRow in range(statesStorage.getSize()):
    #         # Get the current row's StateVector
    #         stateVec = statesStorage.getStateVector(iRow)
            
    #         # Loop through all columns
    #         for iCol in range(columnLabels.getSize()):
    #             colName = columnLabels.get(iCol)
                
    #             # Skip time column
    #             if colName == 'time':
    #                 continue
                
    #             # Check if this is a column that should be excluded from conversion
    #             shouldExclude = False
    #             for excludeCol in excludeColumns:
    #                 if excludeCol in colName:  # Use 'in' to match path-based names or coordinate names
    #                     shouldExclude = True
    #                     break
                
    #             if not shouldExclude:
    #                 # For StateVector, getData() returns the values without the time column
    #                 # So we need to adjust the index (subtract 1 since time is at index 0)
    #                 try:
    #                     # Get current value (note the index adjustment)
    #                     degValue = stateVec.getData().get(iCol - 1)
                        
    #                     # Convert to radians
    #                     radValue = degValue * (np.pi / 180.0)
                        
    #                     # Set the new value
    #                     stateVec.getData().set(iCol - 1, radValue)
                        
    #                 except Exception as e:
    #                     print(f"Error converting column {colName}, index {iCol-1}: {str(e)}")
    #                     continue
    #             else:
    #                 print(f"Excluding {colName} from degree to radian conversion")
        
    #     # Update storage to indicate it's no longer in degrees
    #     statesStorage.setInDegrees(False)
        
    # elif inDegrees and outDegrees:
    #     #Change the storage label back to specifying indegrees=yes
    #     statesStorage.setInDegrees(True)
    # elif not inDegrees and outDegrees:
    #     #Convert radians to degrees
    #     kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesStorage)
    #     #Reset labeling for degrees
    #     statesStorage.setInDegrees(True)
    
    # #Appropriately set output in degrees or radians
    if inDegrees and not outDegrees:
    #     #Convert degrees values to radians for consistency with the current
    #     #file label (defaults back to inDegrees=no). Radians seem to work
    #     #better with the Moco process as well.
    
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertDegreesToRadians(statesTable)
    # elif inDegrees and outDegrees:
    # #     #Change the storage label back to specifying indegrees=yes
    #     statesTable.setInDegrees(True)
    elif not inDegrees and outDegrees:
    #     #Convert radians to degrees
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesTable)
        #Reset labeling for degrees
        try:
            statesTable.addTableMetaData("inDegrees", "yes")
        except:
            pass #in case addTableMetaData isn't in OpenSim version
        # statesTable.addTableMetaData("inDegrees", "yes")
    
    #Write to file
    osim.STOFileAdapter().write(statesTable, outputFileName)


def calculate_plate_torques_from_cop_forces(grf_file_path, output_file_path=None, verbose=True):
    """
    Inputs a grf file with adjusted CoP values to mimic overground walking; am now adjusting torque values to accomodate
    Back-calculates plate torques from forces and COP coordinates in a GRF file.
    
    Coordinate system:
    - X-axis: forward (anterior-posterior)
    - Y-axis: up (vertical)  
    - Z-axis: right (medial-lateral)
    
    Parameters:
    -----------
    grf_file_path : str
        Path to the .sto GRF file containing force and COP data
    output_file_path : str, optional
        Path to save the file with calculated torques. If None, returns data without saving
    verbose : bool, optional
        Print progress information (default: True)
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing original data with recalculated torques
    """
    
    def load_sto_file(filepath):
        """Load OpenSim .sto file into pandas DataFrame"""
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Find header end
        header_end = None
        for i, line in enumerate(lines):
            if 'endheader' in line.lower():
                header_end = i
                break
        
        if header_end is None:
            raise ValueError("Could not find 'endheader' in file")
        
        # Parse headers and data
        headers = lines[header_end + 1].strip().split('\t')
        data_lines = lines[header_end + 2:]
        
        data = []
        for line in data_lines:
            if line.strip():
                data.append([float(x) for x in line.strip().split('\t')])
        
        return pd.DataFrame(data, columns=headers)
    
    def save_sto_file(df, filepath, original_filepath=None):
        """Save DataFrame to OpenSim .sto format"""
        # Read original header if available
        header_lines = []
        if original_filepath:
            with open(original_filepath, 'r') as f:
                lines = f.readlines()
            for line in lines:
                if 'endheader' in line.lower():
                    break
                header_lines.append(line.rstrip())
        else:
            # Create basic header
            header_lines = [
                "version=1",
                f"nRows={len(df)}",
                f"nColumns={len(df.columns)}",
                "inDegrees=no",
                "name=calculated_plate_torques"
            ]
        
        with open(filepath, 'w') as f:
            for line in header_lines:
                f.write(line + '\n')
            f.write('endheader\n')
            f.write('\t'.join(df.columns) + '\n')
            df.to_csv(f, sep='\t', index=False, header=False, float_format='%.6f')
    
    def calculate_plate_torque(cop_x, cop_y, cop_z, force_x, force_y, force_z):
        """
        Calculate plate torques from COP position and ground reaction forces.
        
        For a force plate, the torque represents the moment about the plate origin.
        Using cross product: Torque = COP_vector × Force_vector
        
        Coordinate system: X=forward, Y=up, Z=right
        
        Returns:
        --------
        mx, my, mz : torque components about plate origin
        """
        # Torque = r × F where r is position vector from origin to COP
        mx = cop_y * force_z - cop_z * force_y  # Torque about X-axis (sagittal plane rotation)
        my = cop_z * force_x - cop_x * force_z  # Torque about Y-axis (transverse plane rotation) 
        mz = cop_x * force_y - cop_y * force_x  # Torque about Z-axis (frontal plane rotation)
        
        return mx, my, mz
    
    if verbose:
        print(f"Loading GRF data from: {grf_file_path}")
    
    # Load original GRF data
    grf_data = load_sto_file(grf_file_path)
    
    if verbose:
        print(f"GRF data shape: {grf_data.shape}")
        print(f"Available columns: {list(grf_data.columns)}")
    
    # Create copy for modification
    calculated_data = grf_data.copy()
    
    # Define possible column name patterns for forces, COP, and torques
    # Updated patterns for column format
    force_patterns = [
        ('ground_{}_force_vx', 'ground_{}_force_vy', 'ground_{}_force_vz'),
        ('ground_force_v{}_px', 'ground_force_v{}_py', 'ground_force_v{}_pz'),
        ('ExternalForce_{}_fx', 'ExternalForce_{}_fy', 'ExternalForce_{}_fz'),
        ('force_{}_x', 'force_{}_y', 'force_{}_z'),
        ('f{}_x', 'f{}_y', 'f{}_z')
    ]
    
    cop_patterns = [
        ('ground_{}_force_px', 'ground_{}_force_py', 'ground_{}_force_pz'),
        ('ground_force_v{}_px', 'ground_force_v{}_py', 'ground_force_v{}_pz'),
        ('ExternalForce_{}_px', 'ExternalForce_{}_py', 'ExternalForce_{}_pz'),
        ('cop_{}_x', 'cop_{}_y', 'cop_{}_z'),
        ('px{}', 'py{}', 'pz{}')
    ]
    
    torque_patterns = [
        ('ground_{}_torque_x', 'ground_{}_torque_y', 'ground_{}_torque_z'),
        ('ground_torque_v{}_px', 'ground_torque_v{}_py', 'ground_torque_v{}_pz'),
        ('ExternalForce_{}_mx', 'ExternalForce_{}_my', 'ExternalForce_{}_mz'),
        ('torque_{}_x', 'torque_{}_y', 'torque_{}_z'),
        ('mx{}', 'my{}', 'mz{}')
    ]
    
    # Process each foot/plate (typically numbered 1 and 2)
    for foot_id in [1, 2]:
        side_name = 'right' if foot_id == 1 else 'left'
        
        if verbose:
            print(f"\nProcessing {side_name} foot (ID: {foot_id})...")
        
        # Find force columns
        force_cols = {}
        for pattern in force_patterns:
            fx_col = pattern[0].format(foot_id)
            fy_col = pattern[1].format(foot_id)
            fz_col = pattern[2].format(foot_id)
            
            if all(col in calculated_data.columns for col in [fx_col, fy_col, fz_col]):
                force_cols = {'x': fx_col, 'y': fy_col, 'z': fz_col}
                break
        
        if not force_cols:
            if verbose:
                print(f"  Warning: Could not find force columns for foot {foot_id}")
            continue
        
        # Find COP columns
        cop_cols = {}
        for pattern in cop_patterns:
            px_col = pattern[0].format(foot_id)
            py_col = pattern[1].format(foot_id)
            pz_col = pattern[2].format(foot_id)
            
            if all(col in calculated_data.columns for col in [px_col, py_col, pz_col]):
                cop_cols = {'x': px_col, 'y': py_col, 'z': pz_col}
                break
        
        if not cop_cols:
            if verbose:
                print(f"  Warning: Could not find COP columns for foot {foot_id}")
            continue
        
        # Find torque columns
        torque_cols = {}
        for pattern in torque_patterns:
            mx_col = pattern[0].format(foot_id)
            my_col = pattern[1].format(foot_id)
            mz_col = pattern[2].format(foot_id)
            
            if all(col in calculated_data.columns for col in [mx_col, my_col, mz_col]):
                torque_cols = {'x': mx_col, 'y': my_col, 'z': mz_col}
                break
        
        if not torque_cols:
            if verbose:
                print(f"  Warning: Could not find torque columns for foot {foot_id}")
            continue
        
        # Extract data
        forces = {
            'x': calculated_data[force_cols['x']].values,
            'y': calculated_data[force_cols['y']].values,
            'z': calculated_data[force_cols['z']].values
        }
        
        cop = {
            'x': calculated_data[cop_cols['x']].values,
            'y': calculated_data[cop_cols['y']].values,
            'z': calculated_data[cop_cols['z']].values
        }
        
        if verbose:
            print(f"  Found columns:")
            print(f"    Forces: {list(force_cols.values())}")
            print(f"    COP: {list(cop_cols.values())}")
            print(f"    Torques: {list(torque_cols.values())}")
        
        # Calculate plate torques for each time point
        num_points = len(forces['x'])
        new_torques = {'x': np.zeros(num_points), 'y': np.zeros(num_points), 'z': np.zeros(num_points)}
        
        for i in range(num_points):
            mx, my, mz = calculate_plate_torque(
                cop['x'][i], cop['y'][i], cop['z'][i],
                forces['x'][i], forces['y'][i], forces['z'][i]
            )
            new_torques['x'][i] = mx
            new_torques['y'][i] = my
            new_torques['z'][i] = mz
        
        # Update torque columns in the dataframe
        calculated_data[torque_cols['x']] = new_torques['x']
        calculated_data[torque_cols['y']] = new_torques['y']
        calculated_data[torque_cols['z']] = new_torques['z']
        
        if verbose:
            print(f"  Calculated torques for {num_points} time points")
            print(f"  Torque ranges:")
            print(f"    Mx: [{new_torques['x'].min():.3f}, {new_torques['x'].max():.3f}] N·m")
            print(f"    My: [{new_torques['y'].min():.3f}, {new_torques['y'].max():.3f}] N·m")
            print(f"    Mz: [{new_torques['z'].min():.3f}, {new_torques['z'].max():.3f}] N·m")
    
    # Save calculated data if output path provided
    if output_file_path:
        save_sto_file(calculated_data, output_file_path, grf_file_path)
        if verbose:
            print(f"\nSaved calculated torque data to: {output_file_path}")
    
    if verbose:
        print("\nPlate torque calculation completed successfully")
    
    return calculated_data




# %% ----- End of osimFunctions.py -----


# def load_sto_tsv(filename):
#     """
#     Loads an STO file (tab-delimited) using 'endheader' as a landmark.
    
#     Args:
#         filename (str): Path to the .sto file to load
        
#     Returns:
#         dict or pd.DataFrame: Data with fields named according to headers
        
#     Notes:
#         - The row containing 'endheader' is used as a landmark
#         - The row immediately after 'endheader' contains column headers
#         - The data begins two rows after the 'endheader' row
#     """
#     # Check if file exists (Python will raise FileNotFoundError if it doesn't)
    
#     # Read file line by line to find 'endheader'
#     with open(filename, 'r') as f:
#         lines = f.readlines()
    
#     # Find the line with 'endheader'
#     endheader_idx = None
#     for i, line in enumerate(lines):
#         if 'endheader' in line:
#             endheader_idx = i
#             break
    
#     if endheader_idx is None:
#         raise ValueError(f"No 'endheader' found in file: {filename}")
    
#     # Headers are in the line after 'endheader'
#     header_line = lines[endheader_idx + 1]
    
#     # Data starts two lines after 'endheader'
#     data_start_idx = endheader_idx + 2
    
#     # Parse headers (tab-delimited)
#     headers = [h.strip() for h in header_line.split('\t') if h.strip()]
    
#     # Create a dictionary to store the data
#     data_dict = {header: [] for header in headers}
    
#     # Parse the data
#     for i in range(data_start_idx, len(lines)):
#         line = lines[i].strip()
#         if not line:  # Skip empty lines
#             continue
        
#         values = line.split('\t')
        
#         # Make sure we have the right number of values
#         if len(values) >= len(headers):
#             for j, header in enumerate(headers):
#                 try:
#                     data_dict[header].append(float(values[j]))
#                 except ValueError:
#                     # Handle non-numeric data
#                     data_dict[header].append(values[j])
    
#     # Convert to pandas DataFrame
#     df = pd.DataFrame(data_dict)
    
#     # Display available fields
#     print("Available data fields:")
#     print(df.columns.tolist())
    
#     # Return DataFrame
#     return df




# def process_sto_file(input_file, output_file, inDegrees=None, outDegrees=None):
#     """
#     Process an STO file by loading it, converting angle units if specified,
#     calculating derivatives, and writing to a new file.
    
#     Args:
#         input_file (str): Path to the input .sto file
#         output_file (str): Path to save the processed .sto file
#         inDegrees (bool, optional): If True, input angles are in degrees
#         outDegrees (bool, optional): If True, output angles should be in degrees
    
#     Returns:
#         pd.DataFrame: The processed data
#     """
#     # Load the data using load_sto_tsv
#     df = load_sto_tsv(input_file)
    
#     # Create a copy to work with
#     processed_df = df.copy()
    
#     # Columns that should not be converted
#     non_convertible_columns = ['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz']
    
#     # Get columns to convert (all columns except those in non_convertible_columns)
#     convertible_columns = [col for col in df.columns if col not in non_convertible_columns]
    
#     # Perform angle conversion if requested
#     if inDegrees is True and outDegrees is False:
#         # Convert from degrees to radians
#         for col in convertible_columns:
#             processed_df[col] = df[col] * np.pi / 180.0
#             print(f"Converted {col} from degrees to radians")
#     elif inDegrees is False and outDegrees is True:
#         # Convert from radians to degrees
#         for col in convertible_columns:
#             processed_df[col] = df[col] * 180.0 / np.pi
#             print(f"Converted {col} from radians to degrees")
    
#     # Calculate derivatives (speeds) for all columns except 'time'
#     derivative_columns = [col for col in df.columns if col != 'time']
#     time_series = df['time'].values
    
#     for col in derivative_columns:
#         # Calculate derivative (central difference method)
#         values = df[col].values
        
#         # Calculate derivative using numpy's gradient function
#         # This uses central differences in the interior and one-sided differences at the boundaries
#         derivative = np.gradient(values, time_series)
        
#         # Create new column name for the derivative
#         speed_col = f"jointset/{col}/speed"
        
#         # Add the derivative to the processed dataframe
#         processed_df[speed_col] = derivative
#         print(f"Calculated derivative for {col} and stored as {speed_col}")
    
#     # Rename the original columns to follow the naming convention
#     original_columns = list(df.columns)
#     rename_dict = {}
    
#     for col in original_columns:
#         if col != 'time':
#             rename_dict[col] = f"jointset/{col}/value"
    
#     # Apply the renaming
#     processed_df = processed_df.rename(columns=rename_dict)
    
#     # Write the processed data to a new .sto file
#     write_sto_file(processed_df, output_file, outDegrees)
    
#     print(f"Processed data written to {output_file}")
    
#     return processed_df


# def write_sto_file(df, filename, outDegrees=None):
#     """
#     Write a pandas DataFrame to an STO file format.
    
#     Args:
#         df (pd.DataFrame): The DataFrame to write
#         filename (str): The output file path
#         outDegrees (bool, optional): If True, output angles are in degrees
#     """
#     # Create directory if it doesn't exist
#     os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    
#     # Open the file for writing
#     with open(filename, 'w') as f:
#         # Write header information
#         f.write("OpenSimSTO version=1\n")
#         f.write("name=converted_motion\n")
#         f.write("datacolumns={}\n".format(len(df.columns)))
#         f.write("datarows={}\n".format(len(df)))
#         f.write("range={} {}\n".format(df['time'].min(), df['time'].max()))
        
#         # Add inDegrees metadata based on outDegrees
#         if outDegrees is not None:
#             if outDegrees is True:
#                 f.write("inDegrees=yes\n")
#             else:
#                 f.write("inDegrees=no\n")
        
#         f.write("endheader\n")
        
#         # Write column headers
#         f.write("\t".join(df.columns) + "\n")
        
#         # Write data rows
#         df.to_csv(f, sep='\t', index=False, header=False, float_format='%.8f')

# #To test single file cycle
# statesToKinematics(statesFileName='/home/lisca/biomec/data/20240904_dataset/sub01/ses20240904/sub01_strifs_cond00000_speedp0200_0001_0001_rad.sto', 
#                        osimModelFileName='/home/lisca/biomec/data/20240904_dataset/sub01/ses20240904/sub01_SCALED.osim', 
#                        outputFileName='/home/lisca/biomec/data/20240904_dataset/sub01/ses20240904/sub01_strifs_cond00000_speedp0200_0001_0001_mocokin.sto',
#                        inDegrees = False, outDegrees = True)

if __name__ == "__main__":
    import sys
    import validation as validate
    import retrieve_results as get_results
    
    
    
    input_test = "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/guesses/initial/sub03_strifc_cond00000_speed02000_0001_0001_mi_emg.sto"
    # input_test = "/home/lisca/biomec/data/20250306_dataset/results_pool/sub03_strifc_cond00000_speed02000_0001_0003_ik_move.sto"
    moco_solution_file="/home/lisca/biomec/data/20250306_dataset/results_pool/averages/guesses/sub04_strifc_cond00000_speed02000_20res_02mesh.sto"
    model_file = "example_muscle_adjusted_model.osim" # "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim"
    dont_touch = ['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz']
    converted_file = 'z_example.sto' # "/home/lisca/biomec/data/20250306_dataset/results_pool/sub03_strifc_cond00000_speed02000_0001_0003_rad.sto"
    # kinematicsToStates(kinematicsFileName=input_test, 
    #                    osimModelFileName=model_file, 
    #                    outputFileName=converted_file, 
    #                    inDegrees=True,
    #                    outDegrees=False,
    #                    excludeColumns=dont_touch, 
    #                    include_speeds=False, 
    #                    include_accelerations=False)

    kinematicsToStatesForMoco(mocoSolutionFileName=moco_solution_file, outputFileName='z_kinematics.sto', 
                              inDegrees=False, outDegrees=True, verbose=True)

    # Retrieve simulated GRFs generated by each specific contact sphere
    # path_solution = "/home/lisca/biomec/src/athlete/athlete/moco/example_muscle_driven_2_nospring.sto"
    # path_grf = "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_strifc_cond00000_speed02000_0001_0003_grf_move.xml"
    # path_grf_data = "/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub03_strifc_cond00000_speed02000_grf_move_averages.sto"
    # converted_states = "converted_states.sto" # "/home/lisca/biomec/data/20250306_dataset/results_pool/sub03_strifc_cond00000_speed02000_0001_0003_rad.sto"
    # path_problem = "/home/lisca/biomec/src/athlete/athlete/moco/example_muscle_driven_2_nospring_setup.omoco"
    # solution = osim.MocoTrajectory(path_solution)
    # model = osim.ModelProcessor(model_file).process()
    # model.initSystem()

    # # Calculate plate torques from existing force and COP data
    # calculated_data = calculate_plate_torques_from_cop_forces(
    #     grf_file_path=path_grf_data,
    #     output_file_path=path_grf_data # .replace("sub03_strifc_cond00000_speed02000_grf_move_averages", "averages_grfs_updated")
    # )

    # forceNamesRightFoot = [
    # 'forceset/contactHeel_r',
    # 'forceset/contactLateralRearfoot_r', 
    # 'forceset/contactLateralMidfoot_r',
    # 'forceset/contactLateralToe_r',
    # 'forceset/contactMedialToe_r',
    # 'forceset/contactMedialMidfoot_r'
    # ]

    # # Define contact sphere force names for left foot
    # forceNamesLeftFoot = [
    #     'forceset/contactHeel_l',
    #     'forceset/contactLateralRearfoot_l',
    #     'forceset/contactLateralMidfoot_l', 
    #     'forceset/contactLateralToe_l',
    #     'forceset/contactMedialToe_l',
    #     'forceset/contactMedialMidfoot_l'
    # ]

    # force_paths_ankle_muscles = [
    # '/forceset/gasmed_r',
    # '/forceset/gaslat_r', 
    # '/forceset/soleus_r',
    # '/forceset/tibant_r',
    # '/forceset/gasmed_l',
    # '/forceset/gaslat_l',
    # '/forceset/soleus_l',
    # '/forceset/tibant_l'
    # ]
    
    # # Get net moments from solution
    # solution_trajectory = osim.MocoTrajectory(path_solution)
    # force_paths = osim.StdVectorString()
    # force_paths.push_back(path_grf)
    # # net_moments = calc_net_joint_moments(model, solution_trajectory)
    # # net_moments = study.calcGeneralizedForces(solution_trajectory, force_paths_ankle_muscles)
    # # net_moments_file = path_solution.replace('.sto', '_net_moments.sto')
    # # osim.STOFileAdapter.write(net_moments, net_moments_file)

    # path_comp = path_solution.replace('.sto', 'errors.sto')
    # get_results.analyze_moco_ik_tracking_errors(moco_solution_file=path_solution, ik_kinematics_file=converted_states, output_file=path_comp, 
    #                                 model_file=model_file, omoco_setup_file=path_problem)


    # external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)
    # osim.STOFileAdapter.write(external_loads, os.path.join(path_solution + '_simulated_external_loads.sto'))


    # # path_spring_solution = path_solution.replace('.sto', '_spring_forces_mtp.sto')
    # # get_results.write_mtp_spring_forces(solution=solution, model=model, output_file=path_spring_solution)

    # neg_force_check = path_solution.replace('.sto', '_neg_muscles.csv')
    # validate.get_negative_muscle_forces(solution=solution, model=model, output_path=neg_force_check)

    # # Store the muscle data from the solution
    # muscle_data = validate.calc_muscle_mechanics(model, solution)
    # muscle_data_file = path_solution.replace('.sto', '_muscle_mechanics.csv')
    # osim.STOFileAdapter.write(muscle_data, muscle_data_file)


    # if len(sys.argv) < 3:
    #     print("Usage: python process_sto_file.py <input_file> <output_file> [--deg2rad|--rad2deg]")
    #     sys.exit(1)
    
    # input_file = sys.argv[1]
    # output_file = sys.argv[2]
    
    # inDegrees = None
    # outDegrees = None
    
    # if len(sys.argv) > 3:
    #     if sys.argv[3] == "--deg2rad":
    #         inDegrees = True
    #         outDegrees = False
    #     elif sys.argv[3] == "--rad2deg":
    #         inDegrees = False
    #         outDegrees = True
    
    # process_sto_file(input_file, output_file, inDegrees, outDegrees)


# %%
