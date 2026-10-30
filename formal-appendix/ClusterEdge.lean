/-
SPDX-License-Identifier: Apache-2.0

A uniform numeric edge-width oracle for the reference multiplication
network after clustering. This does not certify the clique-row query.
-/
import NumericGates
import CrossQueryCorrect

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def gateCount (A : ℕ) : ℕ := A + (A + (A + (A + A * A)))

def clusterEdge (A width g u v : ℕ) : ℕ :=
  crossWidth (gateCount A) (gateCluster A g) (gateEdge A width) u v

theorem clusterEdge_polyT {α : Type} [Enc α] {A width g u v : α → ℕ}
    (hA : Unary A) (hw : PolyT width) (hg : PolyT g) (hu : PolyT u) (hv : PolyT v) :
    PolyT fun a => clusterEdge (A a) (width a) (g a) (u a) (v a) := by
  have hN : Unary fun a => gateCount (A a) := by
    exact Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hc : PolyT fun p : α × ℕ => gateCluster (A p.1) (g p.1) p.2 :=
    gateCluster_polyT (PolyT.comp hA.poly PolyT.id'.fst) (PolyT.comp hg PolyT.id'.fst) PolyT.id'.snd
  have he : PolyT fun p : (α × ℕ) × ℕ => gateEdge (A p.1.1) (width p.1.1) p.1.2 p.2 :=
    gateEdge_polyT (PolyT.comp hA.poly PolyT.id'.fst.fst) (PolyT.comp hw PolyT.id'.fst.fst)
      PolyT.id'.fst.snd PolyT.id'.snd
  exact crossWidth_polyT (color := fun a i => gateCluster (A a) (g a) i)
    (edge := fun a i j => gateEdge (A a) (width a) i j) hN hu hv hc he

theorem clusterEdge_eq (s g N k y : ℕ) (hg : 0 < g) (u v : Fin g) :
    clusterEdge (pA s) (fw s) g u.val v.val = (clusterNet s ℤ N k y g hg).c u v := by
  unfold clusterEdge clusterNet
  apply crossWidth_eq_coarsen (mulNet s ℤ N k y)
    (fun a => clMG s g hg ((MG.enum (pA s)).symm a)) (gateCluster (pA s) g) (gateEdge (pA s) (fw s))
  · intro a
    obtain ⟨i, rfl⟩ := (MG.enum (pA s)).surjective a
    simpa using gateCluster_eq s g hg i
  · intro a b
    obtain ⟨i, rfl⟩ := (MG.enum (pA s)).surjective a
    obtain ⟨j, rfl⟩ := (MG.enum (pA s)).surjective b
    simpa [mulNet, GNet.toNet, mulG] using gateEdge_eq s i j

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.clusterEdge_polyT
#print axioms IntegerFactoring.HalfAppendix.clusterEdge_eq
