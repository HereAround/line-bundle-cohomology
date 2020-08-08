#########################################################################
#                                                                       #
# Cohomology of line bundles on toric varieties using the Cech complex. #
#                                                                       #
#########################################################################

from subprocess import Popen, PIPE
import os
import tempfile
import numpy as np

def choose_k_in_s(k, s):
    if k == 0:
        yield []
        return
    
    for i in xrange(0,len(s)-(k-1)):
        for combination in choose_k_in_s(k-1, s[i+1:]):
            yield [s[i]]+combination


def compute_box(rays, divisor):
    dim = len(rays[0])
    # min, max
    boxsize = [(0,0)]*dim
	
    for ray_system in choose_k_in_s(dim, zip(rays, divisor)):
        rays = np.array( [ item[0] for item in ray_system ] )
        a = np.array( [ -item[1] for item in ray_system ] )
        #rays = matrix(CC, [item[0] for item in ray_system])
        #a = matrix(ZZ, [-item[1] for item in ray_system])

        try:
            m = np.linalg.inv( rays )
        except:
        #except ZeroDivisionError:
            # No solution, just try the next collection
            continue

        solution = m * a.transpose()

        for i in range(dim):
            #m_i = int(RR(real(solution[i][0])))
            m_i = int(solution[i][0])
            min_m, max_m = boxsize[i]
            min_m = min(min_m, m_i)
            max_m = max(max_m, m_i)
            boxsize[i] = (min_m, max_m)
		
    for i in range(dim):
        boxsize[i] = (boxsize[i][0]-1,boxsize[i][1]+1)

    return boxsize

def compute_kth_cohomology(rays, cones, divisor, k):
    # Some basic sanity checks.
    assert len(divisor) == len(rays)
    assert all(all(ray < len(rays) for ray in cone) for cone in cones)
    assert len(cones) > len(rays[0])

    # Dimension of the N lattice. We will compute up to C^{dim}, d^{dim}
    dim = len(rays[0])

    # A bit more convenient representation for computing intersections.
    cones = [set(cone) for cone in cones]

    def dot(a,b):
        return sum([i*j for i,j in zip(a,b)])

    # Pass the relevant information to the C code
    file_name = "demofile.txt"
    fd = open(file_name, "w")
    box = compute_box(rays, divisor)
    
    # Definition of the box.
    fd.write(str(len(box))+"\n")
    for interval in box:
        fd.write("%d %d\n" % interval)

    # The rays.
    fd.write(str(len(rays))+"\n")
    for ray in rays:
        for xi in ray:
            fd.write("%d " % (xi,))
        fd.write("\n")

    # The divisor.
    for ai in divisor:
        fd.write("%d " % (ai,))
    fd.write("\n")

    # The cones.
    fd.write(str(len(cones))+"\n")
    for cone in cones:
        fd.write("%d\n" % (len(cone),))
        for ray in cone:
            fd.write("%d " % (ray,))
        fd.write("\n")
    
    fd.close()
    
    cech = Popen(["./cech_cohomology", file_name, str(k)], stdout=PIPE)
    
    output, errors = cech.communicate()
    result = int(output.strip())
    
    os.unlink(fd)
    
    return result


######################################################################
# Example
######################################################################

## P_2
rays = [[1,0], [0,1], [-1,-1]]
cones = [[0,1], [0,2], [1,2]]
D = [1, 0, 0]
print( compute_kth_cohomology(rays,cones,D,0) )
