import numpy as np
def main():
    psi=[]
    with open("test_30_9_25.txt", "r") as f:
        lines =f.readlines()
        atom_positions = lines[1:]
        for atom in atom_positions:
            data =atom.split()
            psi.append(float(data[5]))
            
    average_psi = np.mean(psi)
    print(f"the average psi 6 is {average_psi}")
    
if __name__ == "__main__":
    main()