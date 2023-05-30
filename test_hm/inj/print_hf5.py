import h5py

def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.items():
        print("  {}: {}".format(key, val))

def print_data(name, obj):
    if isinstance(obj, h5py.Dataset):
        print("Dataset '{}':\n{}\n".format(name, obj[()]))

def print_hdf5_data_and_attributes(file_name):
    with h5py.File(file_name, 'r') as f:
        print("Attributes:")
        f.visititems(print_attrs)
        print("\nData:")
        f.visititems(print_data)

if __name__ == "__main__":
    file_name = "m1_480_m2_120_z_7.hdf"
    print_hdf5_data_and_attributes(file_name)
