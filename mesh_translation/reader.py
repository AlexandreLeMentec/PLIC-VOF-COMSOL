import numpy as np

# Function to parse mesh data and compute dr and dz
def parse_comsol_mesh(file_path, decimal_precision=8):
    r_coords = []
    z_coords = []
    mode = 0
    dr = []
    dz = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('# Mesh vertex coordinates'):
                mode = 1
                continue
            elif not line.strip() and mode != 2:
                mode = 0
                continue

            if line.startswith('# Type #2'):
                mode = 2
                continue
            elif mode == 2 and line.startswith('# Elements'):
                mode = 3
                continue

            if line.strip() and not line.startswith('#') and mode == 1:
                try:
                    r, z = map(float, line.strip().split())
                    r_coords.append(r)
                    z_coords.append(z)
                except ValueError:
                    continue
            if mode == 3:
                try:
                    dr_i = 0.
                    dz_i = 0.
                    v1, v2, v3, v4 = map(int, line.strip().split())
                    p1 = (r_coords[v1],z_coords[v1])
                    p2 = (r_coords[v2],z_coords[v2])
                    p3 = (r_coords[v3],z_coords[v3])
                    p4 = (r_coords[v4],z_coords[v4])

                    p = [p1, p2, p3, p4]
                    for i in range(4):
                        for j in range(3):
                            if "%.8f" % p[i][0] == "%.8f" % p[(j+i+1)%4][0] :
                                dz_i = abs(p[i][1] - p[(j+i+1)%4][1])
                            if "%.8f" % p[i][1] == "%.8f" % p[(j+i+1)%4][1] :
                                dr_i = abs(p[i][0] - p[(j+i+1)%4][0])
                    dr.append("%.6f" % dr_i)
                    dz.append("%.6f" % dz_i)
                except ValueError:
                    continue
        return dr,dz

# Example usage
file_path = './mesh_translation/mesh_rec.txt'  # Replace with the actual path
dr, dz = parse_comsol_mesh(file_path)
def format_list(lst, line_length=20):
    return ',\n'.join([','.join(map(str, lst[i:i+line_length])) for i in range(0, len(lst), line_length)])

with open('./mesh_translation/variable_lists.txt', 'w') as file:
    file.write(f'dr_list {{{format_list(dr)}}}\n')
    file.write(f'dz_list {{{format_list(dz)}}}\n')