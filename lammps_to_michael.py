def main():
    coordinates=[]
    with open("dump.LJ", "r") as f:
        lines =f.readlines()
        atom_positions = lines[9:1009]
        for atom in atom_positions:
            data =atom.split()
            x, y = float(data[2]), float(data[3])
            coordinates.append([x,y])

    with open(pos_out, "w",newline="") as csvf:
        writer =csv.writer(csvf)
        writer.writerow()
        for atom_i in coordinates:
            writer.writerow(atom_i,f"{coordinates[atom_i][0]}",f"{coordinates[atom_i][1]}")


if __name__ == "__main__":
    main()