import numpy as np

n = 12

hist = {'000000000000': 3014,
        '000000000001': 75,
        '000000010000': 75,
        '000100000000': 55,
        '000100000001': 3,
        '000100001000': 1,
        '000000010001': 1,
        '000100010000': 42,
        '000100010001': 25,
        '000100010010': 1,
        '000100010011': 2,
        '000100010100': 1,
        '000100010111': 1,
        '000100100100': 1,
        '000100100111': 1,
        '000000010100': 1,
        '000101000010': 1,
        '000101010000': 2,
        '000110010001': 1,
        '000110011001': 1,
        '000110011011': 1,
        '000110011111': 1,
        '000110110111': 1,
        '000111001111': 1,
        '000111010101': 1,
        '000111011011': 3,
        '000111011111': 1,
        '000000000010': 29,
        '000000100000': 13,
        '001000000000': 38,
        '001000000101': 1,
        '001000100000': 2,
        '001000100100': 2,
        '001000100101': 1,
        '001000110100': 2,
        '000000100110': 1,
        '001010000000': 1,
        '001010001000': 1,
        '001011111111': 2,
        '001100000000': 1,
        '001100010001': 1,
        '001100010100': 1,
        '001100010111': 2,
        '001100110101': 4,
        '001100110111': 2,
        '001110011111': 1,
        '001110101111': 1,
        '001110111111': 1,
        '001111011111': 2,
        '001111111011': 1,
        '001111111111': 5,
        '000000000100': 35,
        '000001000000': 41,
        '010000000000': 40,
        '010000000001': 2,
        '010000000100': 1,
        '000001000001': 1,
        '010000010000': 3,
        '000001000100': 1,
        '010001000000': 56,
        '010001000001': 2,
        '010001000100': 2,
        '010001001000': 6,
        '010001001010': 4,
        '010001001011': 1,
        '010001011011': 1,
        '010001011111': 1,
        '000001001000': 2,
        '010010001000': 2,
        '010010001010': 1,
        '010010001101': 1,
        '010011000000': 2,
        '010011000010': 4,
        '010011001000': 26,
        '010011001010': 99,
        '010011001011': 5,
        '010011001100': 2,
        '010011001110': 1,
        '010011010111': 1,
        '010011011010': 4,
        '010011011011': 2,
        '010011101010': 1,
        '000001010000': 1,
        '010100000000': 1,
        '010100010111': 2,
        '010101010000': 1,
        '010101011011': 2,
        '010101011111': 1,
        '010110011000': 1,
        '010110011001': 1,
        '010110011111': 1,
        '010111001010': 2,
        '010111001110': 1,
        '010111001111': 1,
        '010111010110': 1,
        '010111010111': 1,
        '010111011000': 1,
        '010111011010': 2,
        '010111011011': 30,
        '010111011101': 2,
        '010111011110': 2,
        '010111011111': 32,
        '010111101111': 2,
        '010111111101': 1,
        '010111111111': 6,
        '000000000110': 1,
        '011000000000': 1,
        '011001000000': 1,
        '011011001000': 1,
        '011011001010': 3,
        '011011100110': 1,
        '011011101111': 1,
        '011011111111': 4,
        '011101011101': 1,
        '011101100101': 1,
        '011101110111': 2,
        '011101111111': 5,
        '011110101111': 1,
        '011110111111': 5,
        '011111011110': 2,
        '011111011111': 9,
        '011111101011': 1,
        '011111101111': 4,
        '011111110011': 1,
        '011111110111': 3,
        '011111111011': 2,
        '011111111100': 1,
        '011111111101': 4,
        '011111111110': 13,
        '011111111111': 94,
        '000000001000': 16,
        '000010000000': 34,
        '100000000000': 52,
        '100000000001': 1,
        '100000000100': 1,
        '000010000001': 1,
        '100000010000': 2,
        '100000010010': 1,
        '100000100000': 17,
        '100000100010': 1,
        '100000100100': 4,
        '100001000000': 1,
        '000010001000': 10,
        '100010000000': 1,
        '000010001010': 7,
        '100010100000': 1,
        '100011101010': 1,
        '100011101100': 1,
        '000000001001': 1,
        '000010010000': 3,
        '100100000000': 2,
        '100100010000': 1,
        '100100010101': 1,
        '100100100000': 1,
        '100100110101': 2,
        '100100110111': 1,
        '100101110111': 1,
        '100101111111': 1,
        '000010011000': 1,
        '100111110111': 1,
        '100111111101': 1,
        '100111111111': 7,
        '101000000100': 2,
        '101000100000': 38,
        '101000100100': 51,
        '101000100101': 13,
        '101000100110': 15,
        '101000110000': 1,
        '101000110100': 2,
        '101000110101': 1,
        '101000110111': 1,
        '000010101000': 1,
        '101010101110': 1,
        '101010111011': 1,
        '101010111111': 2,
        '101011101101': 1,
        '101011101111': 4,
        '101011111111': 2,
        '101100010101': 5,
        '101100010111': 1,
        '101100100100': 1,
        '101100100101': 5,
        '101100100111': 3,
        '101100110001': 6,
        '101100110011': 1,
        '101100110100': 11,
        '101100110101': 108,
        '101100110110': 2,
        '101100110111': 36,
        '101100111111': 1,
        '101101110101': 1,
        '101101110111': 2,
        '101101111110': 1,
        '101101111111': 3,
        '101110011111': 1,
        '101110110101': 3,
        '101110110111': 1,
        '101110111011': 1,
        '101110111101': 1,
        '101110111110': 7,
        '101110111111': 53,
        '101111011011': 1,
        '101111011111': 2,
        '101111101110': 2,
        '101111101111': 3,
        '101111110011': 1,
        '101111110111': 6,
        '101111111001': 1,
        '101111111011': 7,
        '101111111101': 7,
        '101111111110': 7,
        '101111111111': 122,
        '110000010000': 1,
        '110001001010': 1,
        '000011001000': 3,
        '000011001010': 3,
        '000011001011': 1,
        '110011001000': 1,
        '110011001010': 2,
        '110011011111': 1,
        '110011101110': 1,
        '110011101111': 2,
        '110011111011': 1,
        '110011111111': 4,
        '110101111101': 1,
        '110101111111': 3,
        '110110101011': 1,
        '110110111101': 1,
        '110110111111': 2,
        '110111011111': 4,
        '110111100111': 1,
        '110111101111': 2,
        '110111110111': 2,
        '110111111011': 2,
        '110111111101': 3,
        '110111111110': 9,
        '110111111111': 67,
        '111000000100': 1,
        '111001100110': 1,
        '111001101110': 1,
        '111001101111': 3,
        '111001110000': 1,
        '111001110101': 1,
        '111001110111': 1,
        '111001111011': 1,
        '111001111111': 6,
        '111010111111': 2,
        '111011011110': 1,
        '111011011111': 2,
        '111011100110': 1,
        '111011100111': 1,
        '111011101010': 2,
        '111011101011': 2,
        '111011101100': 1,
        '111011101101': 3,
        '111011101110': 18,
        '111011101111': 42,
        '111011110111': 2,
        '111011111011': 3,
        '111011111100': 1,
        '111011111101': 3,
        '111011111110': 13,
        '111011111111': 81,
        '111100110001': 1,
        '111100110101': 1,
        '111100110111': 1,
        '111100111110': 2,
        '111100111111': 3,
        '111101010110': 1,
        '111101011111': 4,
        '111101100110': 1,
        '111101100111': 1,
        '111101101111': 10,
        '111101110110': 2,
        '111101110111': 15,
        '111101111010': 1,
        '111101111011': 8,
        '111101111101': 1,
        '111101111110': 14,
        '111101111111': 114,
        '111110011111': 1,
        '111110101111': 1,
        '111110110011': 1,
        '111110110111': 1,
        '111110111011': 5,
        '111110111101': 2,
        '111110111110': 8,
        '111110111111': 65,
        '111111001111': 6,
        '111111010110': 1,
        '111111010111': 3,
        '111111011011': 2,
        '111111011101': 3,
        '111111011110': 1,
        '111111011111': 39,
        '111111100101': 1,
        '111111100111': 3,
        '111111101011': 2,
        '111111101101': 5,
        '111111101110': 5,
        '111111101111': 111,
        '111111110011': 2,
        '111111110100': 1,
        '111111110101': 4,
        '111111110110': 5,
        '111111110111': 83,
        '111111111000': 1,
        '111111111001': 2,
        '111111111010': 3,
        '111111111011': 64,
        '111111111100': 11,
        '111111111101': 106,
        '111111111110': 160,
        '111111111111': 2231}

cal_matrices = [np.array([[0.97753906, 0.06420898],
                          [0.02246094, 0.93579102]]),
                np.array([[0.98413086, 0.03210449],
                          [0.01586914, 0.96789551]]),
                np.array([[0.98083496, 0.02575684],
                          [0.01916504, 0.97424316]]),
                np.array([[0.99316406, 0.02490234],
                          [0.00683594, 0.97509766]]),
                np.array([[0.97192383, 0.0456543],
                          [0.02807617, 0.9543457]]),
                np.array([[0.98925781, 0.02246094],
                          [0.01074219, 0.97753906]]),
                np.array([[0.98193359, 0.02893066],
                          [0.01806641, 0.97106934]]),
                np.array([[0.98925781, 0.03625488],
                          [0.01074219, 0.96374512]]),
                np.array([[0.98779297, 0.02539062],
                          [0.01220703, 0.97460938]]),
                np.array([[0.99353027, 0.02600098],
                          [0.00646973, 0.97399902]]),
                np.array([[0.99145508, 0.03491211],
                          [0.00854492, 0.96508789]]),
                np.array([[0.98364258, 0.0345459],
                          [0.01635742, 0.9654541]])]