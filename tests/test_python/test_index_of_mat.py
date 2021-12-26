def index_of_mat(state, pos_clbits):
    sub_state = ""
    for pos in pos_clbits:
        sub_state += state[pos]
    print(sub_state)
    return int(sub_state, 2)

print(index_of_mat("01101000", [1,3]))