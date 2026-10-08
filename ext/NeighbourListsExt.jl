


module NeighbourListsExt

using Unitful
import AtomsBase
using NeighbourLists
import EquivariantTensors as ET
using DecoratedParticles: PState
using StaticArrays

function ET.Atoms.interaction_graph(sys::AtomsBase.AbstractSystem, rcut)
   nlist = NeighbourLists.PairList(sys, rcut, int_type = Int)
   return ET.Atoms.nlist2graph(nlist, sys)
end

function ET.Atoms.nlist2graph(nlist::NeighbourLists.PairList, sys::AtomsBase.AbstractSystem)
   ii = copy(nlist.i)
   jj = copy(nlist.j)
   first = copy(nlist.first)
   R_ij = [ NeighbourLists._getR(nlist, n) for n = 1:length(ii) ]
   S_i = [ AtomsBase.species(sys, i) for i in ii ]
   S_j = [ AtomsBase.species(sys, j) for j in jj ]
   X_ij = [ PState(𝐫 = 𝐫, z0 = si, z1 = sj, 𝐒 = shift)
            for (𝐫, si, sj, shift) in zip(R_ij, S_i, S_j, nlist.S) ]

   # for node data we use _only_ the atomic species for now so that we
   # don't even give the option of using position information directly.
   # ... until we sort out how to best handle this in ET.
   X_i = [ PState(𝐫 = ustrip.(AtomsBase.position(sys, i)),
                  z = AtomsBase.species(sys, i))
           for i = 1:length(sys) ]

   cell_vecs_u = AtomsBase.cell_vectors(sys)
   cell_vecs = ntuple( i -> ustrip.(cell_vecs_u[i]),
                       length(cell_vecs_u) )

   sys_data = ( pbc = AtomsBase.periodicity(sys),
               cell = cell_vecs
              )

   G = ET.ETGraph(ii, jj;
                  edge_data = X_ij,
                  node_data = X_i,
                  graph_data = sys_data)
   @assert G.first == first

   return G
end

function ET.Atoms.forces_from_edge_grads(sys::AtomsBase.AbstractSystem, G::ET.ETGraph, ∇E_edges)

   TFRC = typeof(∇E_edges[1].𝐫)
   F = zeros(TFRC, length(sys))

   for (i, j, e) in zip(G.ii, G.jj, ∇E_edges)
      F[i] -= e.𝐫
      F[j] += e.𝐫
   end

   return F
end

end
