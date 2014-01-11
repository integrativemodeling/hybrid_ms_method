sampling.py
# This programme was implemented within the IMP software package 
# The full programme and documentation as well as instructions how to run # it can be found in: # https://github.com/integrativemodeling/hybrid_ms_method
import IMP
import IMP.atom
import IMP.container
import IMP.display
import IMP.statistics
#import IMP.example
import sys, math, os, optparse

from optparse import OptionParser


# Pdb files ro parse as arguments
if len(sys.argv) !=4:
   sys.exit(" Must provide two pdb files and the number of models you want to generate")

# Convert the arguments into strings and number
Firstpdb = str(sys.argv[1])
Secondpdb = str(sys.argv[2])
models = float(sys.argv[3])

#***************************************** 

# the spring constant to use, it doesnt really matter
k=100
# the target resolution for the representation, this is used to specify how detailed
# the representation used should be
resolution=300
# the box to perform everything 
bb=IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(0,0,0),
                             IMP.algebra.Vector3D(300, 300, 300))


# this function creates the molecular hierarchies for the various involved proteins
def create_representation():
    m= IMP.Model()
    all=IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    all.set_name("the universe")
    # create a protein, represented as a set of connected balls of appropriate
    # radii and number, chose by the resolution parameter and the number of
    # amino acids.
    def create_protein(name, ds):
        h=IMP.atom.create_protein(m, name, resolution, ds)
        leaves= IMP.atom.get_leaves(h)
        # for convenience, have one molecular hierarchy containing all molecules
        all.add_child(h)
        r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(c)\
                                                  for c in h.get_children()],
                                                 k)
        if r:
            m.add_restraint(r)
            # only allow the particles to penetrate or separate by 1 angstrom
            m.set_maximum_score(r, k)

    def create_protein_from_pdbs(name, files):
        def create_from_pdb(file):
            sls=IMP.SetLogState(IMP.NONE)
            #datadir = '/home/elina/IMP/kernel/examples/CL/'
            #datadir = '/Users/apolitis/imp-projects/imp-work/'
            #datadir = '/u/cvr/apolitis/docs-imp/package-python/'
            #file =sys.argv[0]
            datadir = os.getcwd()
            print datadir
	    t=IMP.atom.read_pdb( datadir+'/' + file, m,
                                 IMP.atom.ATOMPDBSelector())
            del sls
            #IMP.atom.show_molecular_hierarchy(t)
            c=IMP.atom.Chain(IMP.atom.get_by_type(t, IMP.atom.CHAIN_TYPE)[0])
            if c.get_number_of_children()==0:
                IMP.atom.show_molecular_hierarchy(t)
            # there is no reason to use all atoms, just approximate the pdb shape instead
            s=IMP.atom.create_simplified_along_backbone(c,
                                                        resolution/300.0)
            IMP.atom.destroy(t)
            # make the simplified structure rigid
            rb=IMP.atom.create_rigid_body(s)
#            rb=IMP.atom.create_rigid_body(c)
            rb.set_coordinates_are_optimized(True)
            return s
#            return c
        if len(files) >1:
            p= IMP.Particle(m)
            h= IMP.atom.Hierarchy.setup_particle(p)
            h.set_name(name)
            for i, f in enumerate(files):
                c=create_from_pdb(f)
                h.add_child(c)
                c.set_name(name+" chain "+str(i))
            r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(c)\
                                                      for c in h.get_children()],
                                                     k)
            if r:
                m.add_restraint(r)
                display_restraints.append(r)
                m.set_maximum_score(r, k)
        else:
            h= create_from_pdb(files[0])
            h.set_name(name)
        all.add_child(h)

    create_protein_from_pdbs("A", [Firstpdb])
    create_protein_from_pdbs("B", [Secondpdb])
    #create_protein_from_pdbs("C", ["rpt3_imp.pdb"])
    return (m, all)



# create the needed restraints and add them to the model
def create_restraints(m, all):
    def add_connectivity_restraint(s):
        tr= IMP.core.TableRefiner()
        rps=[]
        for sc in s:
            ps= sc.get_selected_particles()
            rps.append(ps[0])
            tr.add_particle(ps[0], ps)
        # duplicate the IMP.atom.create_connectivity_restraint functionality
        score= IMP.core.KClosePairsPairScore(IMP.core.HarmonicSphereDistancePairScore(0,1),
                                             tr)
        r= IMP.core.MSConnectivityRestraint(score)
        iA = r.add_type([rps[0]])
        iB = r.add_type([rps[1]])
        #iC = r.add_type([rps[2]])
        n1 = r.add_composite([iA, iB])
        #n2 = r.add_composite([iA, iB], n1)
        #n3 = r.add_composite([iA, iC], n1)


        m.add_restraint(r)
        m.set_maximum_score(r, k)

    def add_distance_restraint(s0, s1):
        r=IMP.atom.create_distance_restraint(s0,s1, 0, k)
        m.add_restraint(r)
        # only allow the particles to separate by one angstrom
        m.set_maximum_score(r, k)

    def add_distance_restraint13(s0, s1):
        d13=20.0
	r=IMP.atom.create_distance_restraint(s0,s1, d13, k)
        m.add_restraint(r)
        # only allow the particles to separate by one angstrom
        m.set_maximum_score(r, 4*k)

    evr=IMP.atom.create_excluded_volume_restraint([all])
    m.add_restraint(evr)
    # a Selection allows for natural specification of what the restraints act on
    S= IMP.atom.Selection
    sA=S(hierarchy=all, molecule="A")
    sB=S(hierarchy=all, molecule="B")
    #sC=S(hierarchy=all, molecule="C")
    add_connectivity_restraint([sA, sB])
    #for l in IMP.atom.get_leaves(all):
    #    r= IMP.example.ExampleRestraint(l, k)
    #    m.add_restraint(r)
    #    m.set_maximum_score(k)
#    add_distance_restraint13(S(hierarchy=all, molecule="A", residue_indexes=[(210, 210)]), S(hierarchy=all, molecule="B", residue_indexes=[(220, 220)]))


# find acceptable conformations of the model
def get_conformations(m):
    sampler= IMP.core.MCCGSampler(m)

    sampler.set_bounding_box(bb)
    # magic numbers, experiment with them and make them large enough for things to work
    sampler.set_number_of_conjugate_gradient_steps(150)
    sampler.set_number_of_monte_carlo_steps(20)
    sampler.set_number_of_attempts(models)
    # We don't care to see the output from the sampler
    #sampler.set_log_level(IMP.SILENT)

    # return the IMP.ConfigurationSet storing all the found configurations that
    # meet the various restraint maximum scores.
    cs= sampler.get_sample()
    return cs

# cluster the conformations and write them to a file
def analyze_conformations(cs, all, gs):
    # we want to cluster the configurations to make them easier to understand
    # in the case, the clustering is pretty meaningless
    embed= IMP.statistics.ConfigurationSetXYZEmbedding(cs,
                 IMP.container.ListSingletonContainer(IMP.atom.get_leaves(all)), True)
    cluster= IMP.statistics.create_lloyds_kmeans(embed, 10, 10000)
    # dump each cluster center to a file so it can be viewed.
    for i in range(cluster.get_number_of_clusters()):
        center= cluster.get_cluster_center(i)
        cs.load_configuration(i)
        w= IMP.display.PymolWriter("cluster.%d.pym"%i)
        for g in gs:
            w.add_geometry(g)



# now do the actual work
(m,all)= create_representation()
IMP.atom.show_molecular_hierarchy(all)
create_restraints(m, all)


# in order to display the results, we need something that maps the particles onto
# geometric objets. The IMP.display.Geometry objects do this mapping.
# IMP.display.XYZRGeometry map an IMP.core.XYZR particle onto a sphere
gs=[]
for i in range(all.get_number_of_children()):
    color= IMP.display.get_display_color(i)
    n= all.get_child(i)
    name= n.get_name()
    g= IMP.atom.HierarchyGeometry(n)
    g.set_color(color)
    gs.append(g)

cs= get_conformations(m)

print "found", cs.get_number_of_configurations(), "solutions"

for i in range(0, cs.get_number_of_configurations()):
        cs.load_configuration(i)
        # print the configuration
        print "solution number: ",i,"scored :", m.evaluate(False)

ListScores = []
for i in range(0, cs.get_number_of_configurations()):
        cs.load_configuration(i)
        # print the configuration
        print "solution number: ",i,"scored :", m.evaluate(False)
        ListScores.append(m.evaluate(False))
        print ListScores

f1 = open("out_scores.csv", "w")
f1.write("\n".join(map(lambda x: str(x), ListScores)))
f1.close()

# for each of the configuration, dump it to a file to view in pymol
for i in range(0, cs.get_number_of_configurations()):
    cs.load_configuration(i)
    w= IMP.display.PymolWriter("config.%d.pym"%i)
    for g in gs:
        w.add_geometry(g)

analyze_conformations(cs, all, gs)





