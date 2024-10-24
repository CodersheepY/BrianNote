import sys
import os
import numpy as np 

print('Please read the head part of this script and get more information!')
print("""
###################################
#                                 #
# for VASP 5.2 or higher versions #
#                                 #
###################################
""")

if len(sys.argv) <= 1:
    print('\n' + ' Warning! ' * 3 + '\n')
    print('You did not select the input file to be converted. \n By default, we are going to convert your CONTCAR.\n')
    
    if not os.path.isfile("POSCAR") and not os.path.isfile("CONTCAR"):
        print("Error:" * 3 + "\n Can not find either POSCAR or CONTCAR!!\n")
        exit()
    else:
        if os.path.isfile("CONTCAR"):
            script = sys.argv[0]
            file_to_be_converted = "CONTCAR"
            print("\n Conversion starts........\n")
        else:
            script = sys.argv[0]
            file_to_be_converted = "POSCAR"
            print("\n There is no CONTCAR in your directory. \n \n POSCAR Conversion starts........\n")
            
if len(sys.argv) == 2:
    print("\n%s Conversion starts......" % sys.argv[1]) 
    script, file_to_be_converted = sys.argv[:2]
else:
    print("\n%s Conversion starts......\n" % sys.argv[1]) 
    script, file_to_be_converted, fixedlayer = sys.argv[:3]
    fixedlayer = int(fixedlayer)

def get_infor():
    with open(file_to_be_converted, 'r') as f:
        lines = f.readlines()
    
    num_atoms = sum([int(x) for x in lines[6].split()])
    
    if lines[7][0].lower() == 's':  # With Selected T T T, coordination starts from line 9
        start_num = 9
        is_direct = lines[8][0].lower() == 'd'
    else: 
        start_num = 8
        print('----------------------------------------------------')
        print(f'Pay Attention! There is no TTT in {file_to_be_converted}')
        print('----------------------------------------------------')
        is_direct = lines[7][0].lower() == 'd'
    
    a, b, c = [], [], []
    if is_direct:         
        for i in range(2, 5):
            line = [float(i) for i in lines[i].strip().split()]
            a.append(line[0])
            b.append(line[1])
            c.append(line[2])
        vector = np.array([a, b, c])
    else: 
        vector = np.array([[1, 0 , 0], [0, 1, 0], [0, 0, 1]])
    
    return vector, lines, start_num, num_atoms, is_direct

def determinelayers(z_cartesian):
    threshold = 0.5 
    seq = sorted(z_cartesian)
    min_value = seq[0]
    layerscount = {}
    sets = [min_value]
    
    for j in range(len(seq)):
        if abs(seq[j] - min_value) >= threshold:
            min_value = seq[j]
            sets.append(min_value)

    for i in range(1, len(sets) + 1):
        layerscount[i] = []            
        for k in range(len(z_cartesian)):   
            if abs(z_cartesian[k] - sets[i-1]) <= threshold:
                layerscount[i].append(k)

    return layerscount

def convert():
    x_cartesian, y_cartesian, z_cartesian, tf = [], [], [], []
    
    for i in range(start_num, num_atoms + start_num):
        line_data = [float(ele) for ele in lines[i].split()[0:3]]
        line_data = np.array([line_data])
        x, y, z = [sum(k) for k in line_data * vector]     
        x_cartesian.append(x)
        y_cartesian.append(y)
        z_cartesian.append(z)

        if start_num == 9:  # If T T T exist, the start_num will be 9        
            tf.append(lines[i].split()[3:])
        else:
            tf.append(' ')   # If there is no T T T, use space instead. 

    layerscount = determinelayers(z_cartesian)

    with open(file_to_be_converted + '_C', 'w') as file_out:
        for i in range(7):
            file_out.write(lines[i].rstrip() + '\n')  # First 7 lines are kept the same
        
        print('\nFound %s layers!' % len(layerscount))

        if len(sys.argv) >= 3:  # This means that the number for fixing layers is given by the user.
            file_out.write('Selective\n')  
            file_out.write('Cartesian\n') 
            for i in range(1, len(layerscount) + 1):
                if i <= fixedlayer: 
                    for j in layerscount[i]:
                        tf[j] = ['F', 'F', 'F']
                else:
                    for k in layerscount[i]:
                        tf[k] = ['T', 'T', 'T']                    
        else:
            if start_num == 9:  # 9 means there are T T T or F F F in the file
                file_out.write('Selective\n')
                file_out.write('Cartesian\n')        
            else:
                file_out.write('Cartesian\n') 

        for i in range(len(x_cartesian)):
            file_out.write(f"\t{float(x_cartesian[i]):+3.10f}   {float(y_cartesian[i]):+3.10f}   {float(z_cartesian[i]):+3.10f}  {' '.join(tf[i])}\n")

vector, lines, start_num, num_atoms, is_direct = get_infor()

if is_direct: 
    print(f"\n{file_to_be_converted} has Direct Coordinates, Conversion starts.... ")
    convert()
else:
    print(f"\n{file_to_be_converted} has Cartesian Coordinates Already! We are going to fix layers only.")
    convert()

print('-----------------------------------------------------\n')
print(f'\n{file_to_be_converted} with Cartesian Coordinates is named as {file_to_be_converted}_C\n')
print('-----------------------------------------------------\n')
