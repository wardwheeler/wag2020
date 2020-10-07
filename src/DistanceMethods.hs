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

-- | neighborJoining takes a list of leaves and a distance matrixx and returns 
-- an NJ tree
neighborJoining :: V.Vector String -> M.Matrix Double -> String -> TreeWithData
neighborJoining leafNames distMatrix outgroup =
	if M.null distMatrix then error "Null matrix in neighborJoining"
	else
		-- get intial matrices
		let initialBigDMatrix = makeInitialDMatrix distMatrix 0 0 []
			(vertexVect, edgeVect, dMatrix) = addTaxa distMatrix initialBigDMatrix leafNames
			newickString = convertToNewick leafNames outgroup (vertexVect, edgeVect)
			treeCost = getTreeCost (vertexVect, edgeVect)
		in
		(newickString, (vertexVect, edgeVect), treeCost, dMatrix)

-- | makeInitialDMatrix makes adjusted matrix (D) from observed (d) values
-- assumes matrix is square and symmetrical
makeInitialDMatrix :: M.Matrix Double -> Int -> Int -> [(Int, Int, Double)] -> M.Matrix Double
makeInitialDMatrix inObsMatrix row column updateList =
    if M.null inObsMatrix then error "Null matrix in makeInitialDMatrix"
    else if row == M.rows inObsMatrix then M.updateMatrix inObsMatrix updateList
    else if column == M.cols inObsMatrix then makeInitialDMatrix inObsMatrix (row + 1) column updateList
    else 
        let dij = inObsMatrix M.! (row, column)
            divisor = (1.0 / ((fromIntegral $ M.rows inObsMatrix) - 2))
            ri  = ((sum $ M.getFullRow inObsMatrix row) - dij) / divisor
            rj  = ((sum $ M.getFullRow inObsMatrix column) - dij) / divisor
            bigDij = dij - (ri + rj)
            newUpdateList = (row, column, bigDij) : updateList
        in
        makeInitialDMatrix inObsMatrix row (column + 1) newUpdateList

-- | pickNearestUpdateMatrix takes d and D matrices, pickes nearest based on D
-- then updates d and D to reflect new node and distances created
-- updates teh column/row for vertices that are joined to be infinity so
-- won't be chosen to join again
pickNearestUpdateMatrix :: M.Matrix Double -> M.Matrix Double -> (M.Matrix Double, M.Matrix Double, Vertex, Edge, Edge)
pickNearestUpdateMatrix littleDMatrix bigDMatrix =
	if M.null littleDMatrix then error "Null d matrix in pickNearestUpdateMatrix"
	else M.null bigDMatrix then error "Null D matrix in pickNearestUpdateMatrix"
	else
		let (iMin, jMin, distIJ) = getMatrixMinPair bigDMatrix (-1, -1, NT.infinity) 0 0
		in
		if distIJ == NT.infinity then error "No minimum found in pickNearestUpdateMatrix"
		else
			let newVertIndex = M.rows littleDMatrix
				dij = littleDMatrix M.! (row, column)
				bottom = (fromIntegral $ M.rows littleDMatrix) - 2
            	divisor = 1.0 / bottom
            	ri  = ((sum $ M.getFullRow littleDMatrix row) - dij) / divisor
            	rj  = ((sum $ M.getFullRow littleDMatrix column) - dij) / divisor
				diMinNewVert = (dij / 2.0) - ((ri - rj) / (2.0 * bottom))
				djMinNewVert = dij - diMinNewVert

				-- get distances to existing vertices
				vertDistList = replicate (M.rows littleDMatrix) newVertIndex
				otherVertList = [0..(M.rows littleDMatrix)]
				pairVertList = (zip vertDistList otherVertList) ++ [(newVertIndex,newVertIndex)]
				newLittleDRow = fmap getNewDist pairVertList
				newLittleDMatrix = M.addMatrixRow littleDMatrix newLittleDMatrix
				newBigDRow = fmap (makeBigDRow newLittleDMatrix) pairVertList

				-- change D valeus to INFTY so won'e be chosen as pairs
				iIFTYList = zip3 (replicate (M.rows newLittleDMatrix) iMin) [0..(M.rows newLittleDMatrix)] (replicate (M.rows newLittleDMatrix) NT.infinity)
				jIFTYList = zip3 (replicate (M.rows newLittleDMatrix) jMin) [0..(M.rows newLittleDMatrix)] (replicate (M.rows newLittleDMatrix) NT.infinity)
				newBigDMatrix = M.updateMatrix (M.addMatrixRow bigDMatrix newBigDRow) (iIFTYList ++ jIFTYList)

				-- create new edges
				newEdgeI = (newVertIndex, iMin, diMinNewVert)
				newEdgeJ = (newVertIndex, jMin,djMinNewVert)
			in
			(newLittleDMatrix, newBigDMatrix, newVertIndex, newEdgeI, newEdgeJ)


