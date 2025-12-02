import tensorflow
import torch

def print_info():
    print(f'tensorflow report: {tensorflow.config.list_physical_devices()}')
    print(f'tensorflow report: {tensorflow.test.is_gpu_available()}')
    print(f'pytorch    report: {torch.cuda.is_available()}')

if __name__ == '__main__':
    print_info()
