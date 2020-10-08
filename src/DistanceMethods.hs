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
import qualified Data.Number.Transfinite           as NT

-- | uPGMA takes a list of leaves and a distance matrixx and returns 
-- an UGPGMA tree
uPGMA :: V.Vector String -> M.Matrix Double -> Int -> TreeWithData
uPGMA leafNames distMatrix outgroup =
  emptyTreeWithData

-- | wPGMA takes a list of leaves and a distance matrixx and returns 
-- an wPGMA tree
wPGMA :: V.Vector String -> M.Matrix Double -> Int -> TreeWithData
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
neighborJoining :: V.Vector String -> M.Matrix Double -> Int -> TreeWithData
neighborJoining leafNames distMatrix outgroup =
    if M.null distMatrix then error "Null matrix in neighborJoining"
    else
        -- get intial matrices
        let initialBigDMatrix = makeDMatrix distMatrix [] 0 0 []
            numLeaves = V.length leafNames
            leafVertexVect = V.fromList $ [0..(numLeaves - 1)]
            nJTree = addTaxaNJ distMatrix initialBigDMatrix numLeaves (leafVertexVect, V.empty) []
            newickString = convertToNewick leafNames outgroup nJTree
            treeCost = getTreeCost nJTree
        in
        -- trace ("NJ tree " ++ show nJTree)
        (newickString, nJTree, treeCost, distMatrix)

-- | sumAvail sums the row only thos values not already added to tree
sumAvail :: [Int] -> Int -> [Double] -> Double
sumAvail vertInList index distList =
  if null distList then 0.0
  else if index `elem` vertInList then sumAvail vertInList (index + 1) (tail distList)
  else 
      let firstDist = head distList
      in
      firstDist + (sumAvail vertInList (index + 1) (tail distList))



-- | makeIDMatrix makes adjusted matrix (D) from observed (d) values
-- assumes matrix is square and symmetrical
-- makes values Infinity if already added
  -- adjust ri and rj to bew based on on values not in termInList
makeDMatrix :: M.Matrix Double -> [Int] -> Int -> Int -> [(Int, Int, Double)] -> M.Matrix Double
makeDMatrix inObsMatrix vertInList row column updateList =
    if M.null inObsMatrix then error "Null matrix in makeInitialDMatrix"
    else if row == M.rows inObsMatrix then M.updateMatrix inObsMatrix updateList
    else if column == M.cols inObsMatrix then makeDMatrix inObsMatrix vertInList (row + 1) 0 updateList
    else if column == row then makeDMatrix inObsMatrix vertInList row (column + 1) ((row, column, 0.0) : updateList)
    else if (column `elem` vertInList) || (row `elem` vertInList) then 
      makeDMatrix inObsMatrix vertInList row (column + 1) ((row, column, NT.infinity) : updateList)
    else
        let dij = inObsMatrix M.! (row, column)
            divisor = ((fromIntegral $ M.rows inObsMatrix) - 2) - (fromIntegral $ length vertInList)
            ri  = (sumAvail vertInList 0 $ M.getFullRow inObsMatrix row) 
            rj  = (sumAvail vertInList 0 $ M.getFullRow inObsMatrix column)
            bigDij = dij - ((ri + rj) / divisor)
        in
        makeDMatrix inObsMatrix vertInList row (column + 1) ((row, column, bigDij) : updateList)

-- | pickNearestUpdateMatrix takes d and D matrices, pickes nearest based on D
-- then updates d and D to reflect new node and distances created
-- updates teh column/row for vertices that are joined to be infinity so
-- won't be chosen to join again
pickNearestUpdateMatrix :: M.Matrix Double -> M.Matrix Double -> [Int] -> (M.Matrix Double, M.Matrix Double, Vertex, Edge, Edge, [Int])
pickNearestUpdateMatrix littleDMatrix bigDMatrix vertInList =
    if M.null littleDMatrix then error "Null d matrix in pickNearestUpdateMatrix"
    else if M.null bigDMatrix then error "Null D matrix in pickNearestUpdateMatrix"
    else
        let (iMin, jMin, distIJ) = getMatrixMinPair bigDMatrix (-1, -1, NT.infinity) 0 0
        in
        -- trace ("First pair " ++ show (iMin, jMin, distIJ)) (
        if distIJ == NT.infinity then error "No minimum found in pickNearestUpdateMatrix"
        else
            let -- new vertex is size of distance matrix (0 indexed)
              newVertIndex = M.rows littleDMatrix
              dij = littleDMatrix M.! (iMin, jMin)
              divisor = (fromIntegral $ M.rows littleDMatrix) - 2 - (fromIntegral $ length vertInList)
              -- only count values of those in
              ri  = (sumAvail vertInList 0 $ M.getFullRow littleDMatrix iMin) 
              rj  = (sumAvail vertInList 0 $ M.getFullRow littleDMatrix jMin)
              -- seem reversed compared to examples 
              -- diMinNewVert = (dij / 2.0) - ((ri - rj) / (2.0 * divisor))
              -- djMinNewVert = dij - diMinNewVert
              djMinNewVert = (dij / 2.0) - ((ri - rj) / (2.0 * divisor))
              diMinNewVert = dij - djMinNewVert

              newVertInList = vertInList ++ [iMin, jMin]

              -- get distances to existing vertices
              otherVertList = [0..((M.rows littleDMatrix) - 1)]
              newLittleDRow = fmap (getNewDist littleDMatrix dij iMin jMin diMinNewVert djMinNewVert) otherVertList 
              newLittleDMatrix = M.addMatrixRow littleDMatrix (V.fromList $ newLittleDRow ++ [0.0])
              -- recalculate whole D matrix since new row affects all the original ones  (except those merged)
              -- included vertex values set to infinity so won't be chosen later
              newBigDMatrix = makeDMatrix newLittleDMatrix newVertInList 0 0 []

              -- create new edges
              newEdgeI = (newVertIndex, iMin, diMinNewVert)
              newEdgeJ = (newVertIndex, jMin, djMinNewVert)
            in
            (newLittleDMatrix, newBigDMatrix, newVertIndex, newEdgeI, newEdgeJ, newVertInList)
            --)

-- | getNewDist get ditance of new vertex to existing vertices 
getNewDist :: M.Matrix Double -> Double-> Int -> Int -> Double -> Double -> Int -> Double
getNewDist littleDMatrix dij iMin jMin diMinNewVert djMinNewVert otherVert =
  if otherVert == iMin then diMinNewVert
  else if otherVert == jMin then djMinNewVert
  else 
    let dik = littleDMatrix M.! (iMin, otherVert)
        djk = littleDMatrix M.! (jMin, otherVert)
    in
    (dik + djk - dij) / 2.0

-- | addTaxaNJ recursively calls pickNearestUpdateMatrix untill all internal nodes are created
-- recursively called until all (n - 2) internal vertices are created.
addTaxaNJ :: M.Matrix Double -> M.Matrix Double -> Int -> Tree -> [Int] -> Tree
addTaxaNJ littleDMatrix bigDMatrix numLeaves inTree vertInList = 
  let (vertexVect, edgeVect) = inTree
  in
  -- completed tree, all inrternal vertieces have been created
  -- unrooted so n + n -2
  -- trace ("Tree vertex number: " ++ show (V.length vertexVect)) (
  let progress = show  ((fromIntegral (100 * ((V.length vertexVect) - numLeaves))/fromIntegral (numLeaves - 2)) :: Double)
  in
  trace (takeWhile (/='.') progress ++ "%") (
  if V.length vertexVect == (2 * numLeaves) - 2 then 
    let (iMin, jMin, _) = getMatrixMinPair bigDMatrix (-1, -1, NT.infinity) 0 0 
        lastEdge = (iMin, jMin, littleDMatrix M.! (iMin, jMin))
    in 
    (vertexVect, edgeVect `V.snoc` lastEdge)
    
  else 
    let (newLittleDMatrix, newBigDMatrix, newVertIndex, newEdgeI, newEdgeJ, newVertInList) = pickNearestUpdateMatrix littleDMatrix bigDMatrix vertInList 
        newVertexVect = vertexVect `V.snoc` newVertIndex
        newEdgeVect = edgeVect V.++ (V.fromList [newEdgeI, newEdgeJ])
    in 
    --trace (M.showMatrixNicely newLittleDMatrix ++ "\n" ++ M.showMatrixNicely bigDMatrix) 
    addTaxaNJ newLittleDMatrix newBigDMatrix numLeaves (newVertexVect, newEdgeVect) newVertInList
    )