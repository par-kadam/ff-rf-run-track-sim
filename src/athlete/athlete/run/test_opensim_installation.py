import opensim

def print_info():
    # Create a dummy model.
    model_test = opensim.Model()

    print(f'OpenSim version: {opensim.GetVersionAndDate()}')

if __name__ == '__main__':
    print_info()
