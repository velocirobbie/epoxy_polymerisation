import json
import numpy as np
from forcefield.opls_reader import OPLS_Reader

class MakeBond(object):
    def __init__(self,sim, change_data, 
            new_connections=False,
            search_radius = 6, mask_radius = 20,
             outputfile=0, networkfile=0, Nbonds=20):

        self.sim = sim
        self.a1  = int(change_data['atoms'].keys()[0])
        self.a2  = int(change_data['atoms'].keys()[1])
        self.change_data = change_data
        self.r2 = search_radius**2 
        self.networkfile = networkfile
        
        self.sim.inv_vdw_defs = {v:k for k,v in self.sim.vdw_defs.items()}
        
        if type(new_connections) == str:
            self.new_connections = json.load(open(new_connections,'r'))
            for key in self.new_connections:
                subdict = self.new_connections[key]
                self.new_connections[key] = {int(k):v for k,v in subdict.items()}
        else:
            self.new_connections = {'angles':{},'dihedrals':{}}
        self.opls = OPLS_Reader('forcefield/oplsaa.prm')

        # calculate NxM distance**2 matrix between types a and b
        all1_indexes, all1_ids, xyz1 = self.find_all_atoms_with_type(self.a1)
        all2_indexes, all2_ids, xyz2 = self.find_all_atoms_with_type(self.a2)
        if len(xyz1) == 0 or len(xyz2) == 0: raise Exception # no possible bonds left to make
        print 'type1: ',len(xyz1),' type2: ',len(xyz2)

        distances = self.calc_distance_matrix(xyz1,xyz2)
        distances_flat = np.ndarray.flatten(distances)

        # map sorted distnaces to atoms that represent that distance
        index_of_sorted_flat_dist = np.argsort(distances_flat)
        row_of_sorted_dists = index_of_sorted_flat_dist / len(xyz2)
        col_of_sorted_dists = index_of_sorted_flat_dist - len(xyz2)*row_of_sorted_dists

        distances_flat.sort()

        dist_index = -1 # initalise for loop, doing a range() was much slower
        searching = True
        count = 0
        coords_of_bonds_made = np.empty((0,3))
        while searching:
            print '===new search==== (',count,')'
            
            # setup mapping from atom ids to array index
            self.ids2index = np.zeros((max(sim.ids+1)),dtype=int)
            for i in range(len(sim.ids)):
                self.ids2index[sim.ids[i]] = i

            # reset mapping from distance matrix to atom index
            for i in range(len(xyz1)):
                all1_indexes[i] = self.ids2index[all1_ids[i]]
            for i in range(len(xyz2)):
                all2_indexes[i] = self.ids2index[all2_ids[i]]

            found = False
            while True:
                dist_index += 1
                dist = distances_flat[dist_index]
                if dist > self.r2: break # all possible bonds are now too far apart
                print 'dist: ',np.sqrt(dist)
                i = np.array( [row_of_sorted_dists[dist_index]] )
                j = np.array( [col_of_sorted_dists[dist_index]] )
                i_distances = self.calc_distance_matrix(xyz1[i], coords_of_bonds_made)
                j_distances = self.calc_distance_matrix(xyz2[j], coords_of_bonds_made)

                if (   np.all(i_distances > mask_radius**2) 
                    or np.all(j_distances > mask_radius**2)):
                    found = True
                    search_dists_from_index = dist_index + 1
                    break
                else:
                    print 'failed, nearby bond: ', np.sqrt(i_distances.min())
    
            if found:
                a,b = all1_indexes[int(i)], all2_indexes[int(j)]
                if (self.sim.atom_labels[a] != self.a1 and
                        self.sim.atom_labels[b] != self.a2): raise Exception

                self.make_bond(self.sim,a,b)

                coords_of_bonds_made = np.vstack((coords_of_bonds_made, xyz1[i], xyz2[j]))
                count += 1
                found = True


            if count >= Nbonds: searching = False
            if not found: searching = False
        
        self.bonds_made = count
        print self.a1, self.a2, 'Bonds made: ', count
        json.dump(self.new_connections,open('new_connections.json','w'))
        if outputfile:
            with open(outputfile,'a') as f:
                f.write(str(self.a1)+' '+ str(self.a2)+'\t Bonds made: '
                        + str(count) + '\n')


    def make_bond(self,sim,a,b):
        
        def remove_hydrogens(a,b):
            ha = self.find_hydrogen_on(sim,a)
            hb = self.find_hydrogen_on(sim,b)
            #print 'removing',ha,hb,self.sim.ids[ha],self.sim.ids[hb]
            self.remove_atoms(sim, hb,ha)
            
            def reorder(i):
                flag = 0
                if ha < i: flag += 1 
                if hb < i: flag += 1
                return i - flag
            a = reorder(a)
            b = reorder(b)
            self.ids2index = np.zeros((max(sim.ids+1)),dtype=int)
            for i in range(len(sim.ids)):
                self.ids2index[sim.ids[i]] = i
            return a,b

        a,b = remove_hydrogens(a,b)
        
        
        change_data = self.change_data
        
        typea = sim.atom_labels[a]
        typeb = sim.atom_labels[b]
        sim.bonds = np.vstack((sim.bonds,[sim.ids[a],sim.ids[b]]))
        print "new bond: index ",a,b,', labels',self.sim.ids[a],self.sim.ids[b],' mols',sim.molecules[a],sim.molecules[b]
        if self.networkfile:
            with open(self.networkfile,'a') as f:
                f.write(str(a)+'\t'+str(b)+'\t'+
                        str(sim.molecules[a])+'\t'+
                        str(sim.molecules[b])+'\t'+
                        str(typea)+'\t'+str(typeb)+'\n')

        def change_atom(i,vdwi,charge):
            try: new_label = sim.inv_vdw_defs[ vdwi ]
            except KeyError: 
                new_label = max(sim.vdw_defs.keys())+1
                sim.vdw_defs[new_label] = vdwi
                sim.inv_vdw_defs[vdwi] = new_label
                #print 'new_label',new_label
            sim.atom_labels[i] = new_label
            sim.charges[i] = charge
            sim.masses[new_label] = self.opls.mass['m'][vdwi-1]
        
        vdwa = change_data['atoms'][str(typea)]['vdw']
        vdwb = change_data['atoms'][str(typeb)]['vdw']
        change_atom(a,vdwa,change_data['atoms'][str(typea)]['charge'])
        change_atom(b,vdwb,change_data['atoms'][str(typeb)]['charge'])
         
        a_opls_type = self.opls.vdw_type['type'][
                      self.opls.vdw_type['vdw'].index(
                      change_data['atoms'][str(typea)]['vdw'])] 
        b_opls_type = self.opls.vdw_type['type'][
                      self.opls.vdw_type['vdw'].index(
                      change_data['atoms'][str(typeb)]['vdw'])] 
        bond = [a_opls_type,b_opls_type]
        found_bond = 0
        for i in range(len(self.opls.bond['a1'])):
            opls_bond = [self.opls.bond['a1'][i],self.opls.bond['a2'][i]]
            if bond == opls_bond or bond == list(reversed(opls_bond)):
                new_coeffs = [ self.opls.bond['k'][i],
                               self.opls.bond['r'][i]]
                for coeff in sim.bond_coeffs:
                    coeff_list= [ sim.bond_coeffs[coeff][1],
                                  sim.bond_coeffs[coeff][2] ]
                    if new_coeffs == coeff_list:
                        sim.bond_labels = np.append(sim.bond_labels,coeff)
                        found_bond += 1
        if found_bond == 0: raise Exception(found_bond)



        aneighbours = self.find_neighbours(a)
        bneighbours = self.find_neighbours(b)        
        neighbours = aneighbours + bneighbours
        
        if change_data['neighbours']:
          for n in neighbours:
            typen = sim.atom_labels[n]
            vdw_n = str(sim.vdw_defs[typen])
            if vdw_n in change_data['neighbours']:
                new_vdw = change_data['neighbours'][vdw_n]['vdw']
                new_charge = change_data['neighbours'][vdw_n]['charge']
                change_atom(n,new_vdw,new_charge)

        self.add_new_angles(a,b,aneighbours,bneighbours)
        self.add_new_dihedrals(a,b,aneighbours,bneighbours)

        self.update_connection_labels(sim, a,b)
        
    def update_connection_labels(self, sim, a, b):
        for thing in ['angle','dihedral']:
            qlist = thing+'s'
            qlabels = thing+'_labels'

            connectionsa  = set(np.where((
                           getattr(self.sim,qlist)== sim.ids[a]) )[0])
            connectionsb = set(np.where((
                           getattr(self.sim,qlist)== sim.ids[b]) )[0])
            connections = connectionsa & connectionsb
            for connection in connections:
                atom_labels = []
                for atom in getattr(self.sim,qlist)[connection]:
                    atom_labels += [self.sim.atom_labels[self.ids2index[atom]]]
                opls_vdws = []
                opls_types = []
                for label in atom_labels:
                    opls_vdws += [self.sim.vdw_defs[label] ]
                for vdw in opls_vdws:
                    opls_types += [ self.opls.vdw_type['type'][
                                    self.opls.vdw_type['vdw'].index(vdw)] ]
                found = 0
                #print '==='
                #print atom_labels
                #print opls_vdws
                #print opls_types
                #print '==='
                sim_qlabels = getattr(self.sim,qlabels)
                for new_thing in self.new_connections[qlist]:
                    types = self.new_connections[qlist][new_thing]
                    if ( opls_types == types or
                         opls_types == list(reversed(types)) ):
                            found += 1
                            label = new_thing
                if found == 0:
                        label = max( sim_qlabels ) + 1
                        print label
                        self.new_connections[qlist][label] = opls_types
                sim_qlabels[connection] = label

    def add_new_angles(self,a,b,aneighbours,bneighbours):
        for i in set(aneighbours) - {b}:
            self.sim.angles = np.vstack((self.sim.angles, [self.sim.ids[i],
                                                           self.sim.ids[a],
                                                           self.sim.ids[b]]))
            self.sim.angle_labels = np.append(self.sim.angle_labels,0)
        for i in set(bneighbours) - {a}:
            self.sim.angles = np.vstack((self.sim.angles, [self.sim.ids[a],
                                                           self.sim.ids[b],
                                                           self.sim.ids[i]]))
            self.sim.angle_labels = np.append(self.sim.angle_labels,0)

    def add_new_dihedrals(self,a,b,aneighbours,bneighbours):
        def add(dihedral):
            self.sim.dihedrals = np.vstack((self.sim.dihedrals,
                                 [self.sim.ids[index] for index in dihedral] ))
            self.sim.dihedral_labels = np.append(self.sim.dihedral_labels, 0)

        aneighbours = set(aneighbours) - {b}
        bneighbours = set(bneighbours) - {a}
        #add torsions with a,b in the centre
        for i in aneighbours:
            for j in bneighbours:
                add([i,a,b,j])
        #add torsions with b at one end
        for i in aneighbours:
            ineighbours = set( self.find_neighbours(i) ) - {a}
            for j in ineighbours:
                add([j,i,a,b])
        #add torsions with a at one end
        for i in bneighbours:
            ineighbours = set( self.find_neighbours(i) ) - {b}
            for j in ineighbours:
                add([j,i,b,a])

    def find_hydrogen_on(self, sim, centre):
        neighbours = self.find_neighbours(centre)
        found = False
        for neighbour in neighbours:
            #print 'H search',centre,':',neighbour,sim.masses[sim.atom_labels[neighbour]],sim.atom_labels[neighbour]
            if sim.masses[sim.atom_labels[neighbour]] == 1.008:
                found = True
                break
        if not found: raise Exception(centre, neighbours)
        return neighbour

    def remove_atoms(self, sim, a, b):
        sim.coords = np.delete(sim.coords,[a,b],0)
        
        #Remove connectivity
        connectivity_quantities = ['bond','angle','dihedral','improper']
        for q in connectivity_quantities:
          for i in [a,b]:
            qlist = q+'s'
            qlabels = q+'_labels'
            connections = np.where( getattr(sim, qlist)==sim.ids[i])
            new_array = np.delete( getattr(sim, qlist),connections[0], 0)
            setattr(sim, qlist, new_array)
            new_array = np.delete( getattr(sim, qlabels), connections[0])
            setattr(sim, qlabels, new_array)

        atom_quantities = ['molecules','atom_labels','charges','ids']
        for q in atom_quantities:
            new_array = np.delete( getattr(sim, q),[a,b])
            setattr( sim, q, new_array)
    
    def find_neighbours(self, centre):
        bonds = np.where(self.sim.bonds==self.sim.ids[centre])
        neighbours = []
        for bond in np.transpose(bonds):
            neighbours += [self.ids2index[self.sim.bonds[bond[0],bond[1]-1]] ]
        return neighbours
    
    def calc_distance_matrix(self,xyz1,xyz2):
        bx = self.sim.xhi - self.sim.xlo
        by = self.sim.yhi - self.sim.ylo
        bz = self.sim.zhi - self.sim.zlo
        def check_periodic(vec, length):
            half_length = length / 2
            over = np.array(vec > half_length)
            under= np.array(vec < -half_length)
            return vec + (under * length) - (over * length)
            
        x = xyz2[:,0] - xyz1[:,0,np.newaxis]
        x = check_periodic(x, bx)
        y = xyz2[:,1] - xyz1[:,1,np.newaxis]
        y = check_periodic(y, by)
        z = xyz2[:,2] - xyz1[:,2,np.newaxis]
        z = check_periodic(z, bz)
        return x*x + y*y + z*z

    def find_all_atoms_with_type(self, label):
        atom_indexes = []
        atom_ids = []
        xyz   = []
        for i in range(len(self.sim.atom_labels)):
            if self.sim.atom_labels[i] == label:
                atom_indexes += [i]
                atom_ids += [self.sim.ids[i]]
                xyz   += [self.sim.coords[i]]
        return np.array(atom_indexes), np.array(atom_ids), np.array(xyz)

