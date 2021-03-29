# Takes a numpy.ndarray and prints parameters for debug
def array_inspect(arr, title_string, print_all):
    print('--------------------', title_string, '--------------------')
    print('arr.shape', arr.shape)
    print('type(arr)', type(arr))
    if(print_all):
        print('arr= ', arr)
    print('----------------------------------------------------------')