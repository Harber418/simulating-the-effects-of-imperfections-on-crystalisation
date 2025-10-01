
import csv
def main():
    coordinates=[]
    with open("dump.LJ", "r") as f:
        lines =f.readlines()
        for i in range(1,102):
            index1 = i*1009 +9 
            index2 = index1 + 1000
            set1 = []
            atom_positions = lines[index1:index2]
            for atom in atom_positions:
                data =atom.split()
                x, y = float(data[2]), float(data[3])
                set1.append([x,y])
            coordinates.append(set1)
            

    with open("pos_file", "w",newline="") as csvf:
        writer =csv.writer(csvf, delimiter=" ")
        final_snapshot = coordinates[-1]
        for atom_i in final_snapshot:
            writer.writerow([f"{atom_i[0]}",f"{atom_i[1]}"])


if __name__ == "__main__":
    main()