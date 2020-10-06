{- |
Module      :  DistanceMethods.hs
Description :  Module to calculate distance tree construction methods Neightbor-Joining, UPGMA, and WPGMA
              -- but with added refinement based on 4-point metric
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

I
-}

module DistanceMethods (neighborJoining, uPGMA, wPGMA, doWagnerS, performWagnerRefinement) where

import qualified Data.Vector                       as V 
import           Debug.Trace
import qualified SymMatrix                         as M
import           Types
import           Utilities
import qualified Wagner                            as W

-- | neighborJoining takes a list of leaves and a distance matrixx and returns 
-- an NJ tree
neighborJoining :: V.Vector String -> M.Matrix Double -> String -> TreeWithData
neighborJoining leafNames distMatrix outgroup =
  emptyTreeWithData

-- | uPGMA takes a list of leaves and a distance matrixx and returns 
-- an UGPGMA tree
uPGMA :: V.Vector String -> M.Matrix Double -> String -> TreeWithData
uPGMA leafNames distMatrix outgroup =
  emptyTreeWithData

-- | wPGMA takes a list of leaves and a distance matrixx and returns 
-- an wPGMA tree
wPGMA :: V.Vector String -> M.Matrix Double -> String -> TreeWithData
wPGMA leafNames distMatrix outgroup =
  emptyTreeWithData

-- | pulls dWagner function from module Wagner
doWagnerS :: V.Vector String -> M.Matrix Double -> String -> Int -> String -> [V.Vector Int]-> [TreeWithData]
doWagnerS leafNames distMatrix firstPairMethod outgroup addSequence replicateSequences = 
  W.doWagnerS leafNames distMatrix firstPairMethod outgroup addSequence replicateSequences 

-- | pulls Wagner refinement from Wagner module
performWagnerRefinement :: String -> String -> String -> V.Vector String -> Int -> TreeWithData -> [TreeWithData]
performWagnerRefinement refinement saveMethod keepMethod leafNames outGroup inTree = 
  W.performRefinement refinement saveMethod keepMethod leafNames outGroup inTree