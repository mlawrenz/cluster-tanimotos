from msmbuilder import metrics
import sys
import glob
import os
import sys
import optparse
import pickle
import numpy
from msmbuilder import Trajectory

def assign(matrix, database, cutoff):
    distances=-1*numpy.ones(len(database))
    assignments=-1*numpy.ones(len(database), dtype=int)
    for j in xrange(matrix.shape[0]):
        d=matrix[j,:]
        ind=numpy.argmin(d)
        if not d[ind] < cutoff:
            pass
        else:
            assignments[j] = int(numpy.argmin(d))
            distances[j] = d[assignments[j]]
    return assignments, distances
 
def cluster(matrix, distance_cutoff, cluster_cutoff=None):
    if cluster_cutoff is None and distance_cutoff is None:
        raise ValueError("I need some cutoff criterion! both k and distance_cutoff can't both be none")
    if cluster_cutoff is None and distance_cutoff <= 0:
        raise ValueError("With k=None you need to supply a legit distance_cutoff")
    if distance_cutoff is None:
        # set it below anything that can ever be reached
        distance_cutoff = -1
    if cluster_cutoff is None:
        # set k to be the highest 32bit integer
        cluster_cutoff = sys.maxint

    distance_list = numpy.inf * numpy.ones(matrix.shape[0], dtype=numpy.float32)
    assignments = -1 * numpy.ones(matrix.shape[0], dtype=int)
    distance_cutoff=float(distance_cutoff)

    seed=0
    generator_indices = []

    for i in xrange(cluster_cutoff):
        new_ind = seed if i == 0 else numpy.argmax(distance_list)
        #print "K-centers: Finding generator %i. Will finish when % .4f drops below % .4f" % (i, float(distance_list[new_ind]), float(distance_cutoff))
        print "K-centers: Finding generator %i. Will finish when Tanimoto % .4f is above % .4f" % (i, 2-float(distance_list[new_ind]), 2-float(distance_cutoff))
        if distance_list[new_ind] < distance_cutoff:
            break
        new_distance_list = matrix[new_ind, :]
        updated_indices = numpy.where(new_distance_list < distance_list)[0]
        if updated_indices.size==0:
            break
        distance_list[updated_indices] = new_distance_list[updated_indices]
        assignments[updated_indices] = new_ind
        generator_indices.append(new_ind)
    return numpy.array(generator_indices), numpy.array(assignments), numpy.array(distance_list)

                                                           
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-m', '--matrix', dest='matrix',
                      help='NxN matrix of scores')
    parser.add_option('-d', '--dbase', dest='dbase',
                      help='N-sized database filenames, where indices map to the matrix')
    parser.add_option('-c', '--cutoff', dest='cutoff',
                      help='tanimoto cutoff')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    cutoff=float(options.cutoff)
    dir=os.path.dirname('%s' % options.matrix)
    if dir=="":
        dir="./"
    matrix=numpy.loadtxt('%s' % options.matrix)
    # Convert Tanimoto value to Distance from the Max (2.00)
    matrix=2.0-matrix
    #matrix=clean_matrix(matrix)
    dbase=numpy.loadtxt('%s' % options.dbase, dtype=str)

    print "Clustering scores based on Max Tanimoto = 2.00"
    gens, assignments, distances=cluster(matrix, distance_cutoff=2-cutoff )
    frames=numpy.where(assignments!=-1)[0]
    # convert back to Tanimoto scores, instead of distance from the max=2
    distances[frames]=[round((2-i),2) for i in distances[frames]]
    numpy.savetxt('%s/Assignments-%s.dat' % (dir, cutoff), assignments, fmt='%i')
    numpy.savetxt('%s/Distances-c%s.dat' % (dir, cutoff), distances)
    print "Assigned %s ligands to %s clusters" % (len(frames), len(gens))

    for i in gens:
        frames=numpy.where(assignments==i)[0]

        ohandle=open('%s/Gen%s-Members-c%s.dat' % (dir, i, cutoff), 'w')
        for (name, val) in zip(dbase[frames], distances[frames]):
            ohandle.write('%s\t%s\n' % (name, val))
    numpy.savetxt('%s/Gens-c%s.dat' % (dir, cutoff), gens, fmt='%i')
    numpy.savetxt('%s/Gens-Names-c%s.dat' % (dir, cutoff), dbase[gens], fmt='%s')
