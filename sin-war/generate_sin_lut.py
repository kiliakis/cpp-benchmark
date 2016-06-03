
import numpy as np
import sys
import os

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "You must specify the lookup table size and the output file directory"
        exit(-1)
    size = int(sys.argv[1])
    file_dir = sys.argv[2]

    file = open(os.path.join(file_dir, 'sin_lut.h'), 'w')

    file.write('/*\n')
    file.write(' * sin_lut.h\n')
    file.write(' * Script generated file. DO NOT EDIT.\n')
    file.write(' * This file is used to store a pre calculated sin\n')
    file.write(' * lookup table. User can switch on-off this option in the\n')
    file.write(' * configuration.h file. Using a lookup table implementation\n')
    file.write(' * of sin function can result in performane improvements but\n')
    file.write(' * will definetely limit its precision. The precision is defined\n')
    file.write(' * as PI / (2 * table_size)\n')
    file.write(' */\n')
    file.write('\n\n')
    file.write('#ifndef SIN_LUT_H_\n')
    file.write('#define SIN_LUT_H_\n')
    file.write('\n\n')


    step = np.pi / (2*size)

    file.write('static const unsigned lut_size = ' + str(size+1) + ';\n')
    file.write('static const float lut_step = ' + str(step) + ';\n')
    file.write('\n')

    file.write('static const float sin_lut[lut_size] = {\n\n\t\t')

    for i in range(0, size + 1, 4):
        for j in range(0, 4):
            num = np.sin((i+j) * step)
            if(i+j == size):
                file.write(str(num))
                break
            else:
                file.write(str(num) + ',')
        file.write('\n\t\t')

    file.write('\n};\n')


    file.write('#endif /* SIN_LUT_H_ */\n')
