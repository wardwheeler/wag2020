{- |
Module      :  wag2020.hs 
Description :  Progam to calcualte Wagner distance trees ala Farris 1972 
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

Initial implementation stright Farris 1972
  O(n^3)
  closest/furthest 2 taxa to start
  closest join for addition point
  closest taxon to be added 

To do:
  Create compact Newick or other rep for equality conparison of trees
  
  Newick smaller subtree left, bigger right for figtree output
  
  Test keep options and swapping on multiple input trees from build

  precision for showing branch lengths could be an input option ie 0 for not showing them at all.

  Clean up code--remove n^2 bad stuff
    TBR time complexity seems pretty awful--if corrrect

  confirm time complexity for OTU, SPR, TBR

 
-}

module Main where

import System.IO
import System.Environment
import Debug.Trace
import Data.List
import Data.Maybe
import qualified Data.Vector as V hiding (replicateM)
import Text.ParserCombinators.Parsec
import Data.CSV
import Control.Parallel.Strategies
import Control.Concurrent
-- import Control.Monad (replicateM)
import qualified Control.Monad.Parallel as CMP
import Immutable.Shuffle
import qualified Data.Set as Set
import qualified Data.Graph.Inductive.Graph as G
import qualified Data.GraphViz as GV
import qualified Data.Graph.Inductive.PatriciaTree as P
import Data.GraphViz.Printing
import qualified Data.Text.Lazy as T
import qualified Data.Number.Transfinite as NT
import qualified System.Random as Rand
import qualified System.Random.Shuffle as RandS
import System.IO.Unsafe
-- From matrices package (so 0 based)
-- import Data.Matrix as M

-- Local lower diag, boxed, immutable matrices
import qualified SymMatrix as M


type Vertex = Int
type Weight = Double
type Edge = (Vertex, Vertex, Weight)
type Tree = (V.Vector Vertex,V.Vector Edge)

type TreeWithData = (String, Tree, Double, M.Matrix Double)
type SplitTreeData = (V.Vector Edge,V.Vector Edge, Double, V.Vector Edge, M.Matrix Double)

-- | used for comparing tree costs that are Double
epsilon :: Double
epsilon = 0.000000000000001

-- | precision for branch lengths display
precision :: Int
precision = 8

-- | localRoundtakes a double multiplies by 10^precisoin, rounds to integer then divides
-- by precision
localRound :: Double -> Int -> Double
localRound val places =
  let factor = 10.0 ^^ places
      newVal = val * factor
      roundedVal = round newVal :: Int
  in
  fromIntegral roundedVal / factor

-- | showDouble integer number of string postion )including sign and decimal) and returns string
-- of that length--basically truncation
showDouble :: Int -> Double -> String
showDouble places val = show $ localRound val places -- reverse $ dropWhile (== '0') $ reverse $ take places $ show val

-- | test for equality with epsilon
withinEpsilon :: Double -> Double -> Bool
withinEpsilon a b = a == b
  {-
  if abs (a - b) < epsilon then True
  else False
  -}

-- | functions for triples
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

snd3 :: (a,b,c) -> b
snd3 (_,e,_) = e

thd3 :: (a,b,c) -> c
thd3 (_,_,f) = f

fst4 :: (a,b,c,d) -> a
fst4 (e,_,_,_) = e

snd4 :: (a,b,c,d) -> b
snd4 (_,e,_,_) = e

thd4 :: (a,b,c,d) -> c
thd4 (_,_,e,_) = e

-- | getTreeCost takes Tree and returns cost based on sum of edge weights
getTreeCost :: Tree -> Double
getTreeCost inTree =
  V.sum $ V.map getEdgeCost $ snd inTree

-- | orderEdge takes an Edge and puts high index first then lower
orderEdge :: Edge -> Edge
orderEdge (a,b,w) =
  if a > b then (a,b,w)
  else (b,a,w)

-- | orderTree puts Tree edges in order based on edges
orderTree :: Tree -> Tree
orderTree (leaves, edges) =
  let edgeList = sort $ V.toList edges
  in
  (leaves, V.fromList edgeList)

-- | 
-- Map a function over a traversable structure in parallel
-- Preferred over parMap which is limited to lists
-- Add chunking (with arguement) (via chunkList) "fmap blah blah `using` parListChunk chunkSize rseq/rpar"
-- but would have to do one for lists (with Chunk) and one for vectors  (splitAt recusively)
parmap :: Traversable t => Strategy b -> (a->b) -> t a -> t b
parmap strat f = withStrategy (parTraversable strat).fmap f

-- | seqParMap takes strategy,  if numThread == 1 retuns fmap otherwise parmap and 
seqParMap :: Traversable t => Strategy b -> (a -> b) -> t a -> t b
seqParMap strat f =
  if getNumThreads > 1 then parmap strat f
  else fmap f



myStrategy :: (NFData b) => Strategy b
myStrategy = rdeepseq

-- | getMatrixMinPair takes distMatrix initla pinteger pair and value
-- traverses teh matrix and return minimum distance and index pair
-- if tie takes first
getMatrixMinPair :: M.Matrix Double ->  (Int, Int, Double) -> Int -> Int -> (Int, Int, Double)
getMatrixMinPair distMatrix curBest curRow curColumn
  | curRow == M.rows distMatrix = curBest
  | curColumn == M.cols distMatrix = getMatrixMinPair distMatrix curBest (curRow + 1) 0
  | curColumn == curRow = getMatrixMinPair distMatrix curBest curRow (curColumn + 1)
  | otherwise =
  let (_, _, currentBestDistance) = curBest
  in
  if  distMatrix M.! (curRow, curColumn) < currentBestDistance then
    getMatrixMinPair distMatrix (curRow, curColumn, distMatrix M.! (curRow, curColumn)) curRow (curColumn + 1)
  else getMatrixMinPair distMatrix curBest curRow (curColumn + 1)


-- | getMatrixMaxPair takes distMatrix initla pinteger pair and value
-- traverses teh matrix and return maximum distance and index pair
-- if tie takes first
getMatrixMaxPair :: M.Matrix Double ->  (Int, Int, Double) -> Int -> Int -> (Int, Int, Double)
getMatrixMaxPair distMatrix curBest curRow curColumn
  | curRow == M.rows distMatrix = curBest
  | curColumn == M.cols distMatrix = getMatrixMaxPair distMatrix curBest (curRow + 1) 0
  | curColumn == curRow = getMatrixMaxPair distMatrix curBest curRow (curColumn + 1)
  | otherwise =
  let (_, _, currentBestDistance) = curBest
  in
  if  distMatrix M.! (curRow, curColumn) > currentBestDistance then
    getMatrixMaxPair distMatrix (curRow, curColumn, distMatrix M.! (curRow, curColumn)) curRow (curColumn + 1)
  else getMatrixMaxPair distMatrix curBest curRow (curColumn + 1)



-- | getStartingPair returns starying pair for Wagner build
--  closts mnimal cost pair
--  furthest maximal cost pair
--  random chooses uniformly at random from leaf set
getStartingPair :: String -> M.Matrix Double -> Edge
getStartingPair choiceOpt distMatrix
  | choiceOpt == "closest" = getMatrixMinPair distMatrix (-1, -1, NT.infinity) 0 0
  | choiceOpt == "furthest" = getMatrixMaxPair distMatrix (-1 , -1, 0 :: Double) 0 0
  | choiceOpt == "random" = error "Initial pair option 'random' not yet implemented"
  | otherwise = error ("Initial pair option " ++ choiceOpt ++ " unrecognized.  Must be 'closest', 'furthest', or 'random'")

-- | getBestEdgeTree take list of edge tuples and return trhe one with best addition cost
getBestEdgeTree :: V.Vector (Double, Tree, M.Matrix Double) -> Double -> (Double, Tree, M.Matrix Double) -> (Double, Tree, M.Matrix Double)
getBestEdgeTree edgeTreeList curBestCost curBestResult =
  if V.null edgeTreeList then curBestResult
  else
    let (firstAddCost, _, _) = V.head edgeTreeList
    in
    if firstAddCost < curBestCost then getBestEdgeTree (V.tail edgeTreeList) firstAddCost (V.head edgeTreeList)
    else getBestEdgeTree (V.tail edgeTreeList) curBestCost curBestResult

-- | add leaf to an edge creating new tree with distances and add cost, also augmented distance matrix
-- but this for swap so returns entire new3-edge cost  so not Farris triangle it is sum of three diveded by 2
addToEdgeSwap :: M.Matrix Double -> Int -> Tree -> Int -> Edge -> (Double, Tree, M.Matrix Double)
addToEdgeSwap distMatrix leaf initialTree newLeafIndex inEdge =
  let (eVertex, uVertex, inWeight) = inEdge
      (initialVertexVect, initialEdgeVect) = initialTree
      addCost = ((distMatrix M.! (leaf, eVertex)) + (distMatrix M.! (leaf, uVertex)) - (distMatrix M.! (eVertex, uVertex))) / 2.0
      eVertLeafDist = (distMatrix M.! (leaf, eVertex)) - addCost
      uVertLeafDist = (distMatrix M.! (leaf, uVertex)) - addCost
      newVertexVect = V.snoc initialVertexVect leaf
      newEdges = V.fromList [(leaf,newLeafIndex, addCost),(eVertex, newLeafIndex, eVertLeafDist),(uVertex, newLeafIndex, uVertLeafDist)]
      cleanupEdges = V.filter (/= inEdge) initialEdgeVect
      newEdgeVect = cleanupEdges V.++ newEdges
      newTree = (newVertexVect, newEdgeVect)
      -- add new costs from added vertex to each reamaining leaf
      augmentedDistMatrix = getNewDistMatrix distMatrix addCost eVertLeafDist uVertLeafDist eVertex uVertex leaf
  in
  (addCost + eVertLeafDist + uVertLeafDist - inWeight , newTree, augmentedDistMatrix)


-- | add leaf to an edge creating new tree with distances and add cost, also augmented distance matrix
addToEdge :: M.Matrix Double -> Int -> Tree -> Int -> Edge -> (Double, Tree, M.Matrix Double)
addToEdge distMatrix leaf initialTree newLeafIndex inEdge =
  -- trace ("In addToEdge with " ++ (show (leaf, initialTree, newLeafIndex, (M.rows distMatrix), inEdge))) (
  let (eVertex, uVertex, _) = inEdge
      (initialVertexVect, initialEdgeVect) = initialTree
      addCost = ((distMatrix M.! (leaf, eVertex)) + (distMatrix M.! (leaf, uVertex)) - (distMatrix M.! (eVertex, uVertex))) / 2.0
      eVertLeafDist = (distMatrix M.! (leaf, eVertex)) - addCost
      uVertLeafDist = (distMatrix M.! (leaf, uVertex)) - addCost
      newVertexVect = V.snoc initialVertexVect leaf
      newEdges = V.fromList [(leaf,newLeafIndex, addCost),(eVertex, newLeafIndex, eVertLeafDist),(uVertex, newLeafIndex, uVertLeafDist)]
      cleanupEdges = V.filter (/= inEdge) initialEdgeVect
      newEdgeVect = V.map orderEdge $ cleanupEdges V.++ newEdges
      newTree = (newVertexVect, newEdgeVect)
      -- add new costs from added vertex to each reamaining leaf
      augmentedDistMatrix = getNewDistMatrix distMatrix addCost eVertLeafDist uVertLeafDist eVertex uVertex leaf
  in
  (addCost, newTree, augmentedDistMatrix)


-- | addTaxonToTree takes distMatrix, an initialTree, Vector of leavesToAdd, and leaf index to add
-- and retursn a tuple wiht the addition cost, the new tree, the new leaves to add, and new distance matrix (enhanced)
addTaxonToTree :: M.Matrix Double -> Tree -> V.Vector Int -> Int -> Int -> (Double, Tree, V.Vector Int, M.Matrix Double)
addTaxonToTree distMatrix initialTree leavesToAdd newVertexIndex leaf =
  if V.null leavesToAdd then (0.0, initialTree, leavesToAdd, distMatrix)
  else
    let leavesRemaining = V.filter (/= leaf) leavesToAdd
        (_, edgesInitial) = initialTree

        -- Parallelize heretoo much and destroys lazy matrix update
        addEdgeList = V.map (addToEdge distMatrix leaf initialTree newVertexIndex) edgesInitial
        (firstAddCost, _, _) = V.head addEdgeList -- this to initialize getBestEdge below
    in
    --filter for best addition point
    let (addCost, newTree, augmentedDistMatrix) = getBestEdgeTree (V.tail addEdgeList) firstAddCost (V.head addEdgeList)
    in
    (addCost, newTree, leavesRemaining, augmentedDistMatrix)

-- | getBestLeafAdd chooses best leaf to add based on cost field
getBestLeafAdd ::  V.Vector (Double, Tree, V.Vector Int, M.Matrix Double) -> Double -> (Double, Tree, V.Vector Int, M.Matrix Double) -> (Double, Tree, V.Vector Int, M.Matrix Double)
getBestLeafAdd addPosVect curBestCost curBestLeaf =
  if V.null addPosVect then curBestLeaf
  else
    let (thisCost, _, _, _) = V.head addPosVect
    in
    if thisCost < curBestCost then getBestLeafAdd (V.tail addPosVect) thisCost (V.head addPosVect)
    else getBestLeafAdd (V.tail addPosVect) curBestCost curBestLeaf


-- | wagBest takes distMatrix, and intial tree of two leaves, a vector of leavesToAdd, the input nuber of leaves
-- and returns the Farris 1972 distance Wagner adding the "closest" leaf at each iteration
wagBest :: M.Matrix Double -> Tree -> V.Vector Int -> Int -> Int ->  V.Vector Int -> String -> (Tree, V.Vector Int, M.Matrix Double)
wagBest distMatrix inTree leavesToAdd nOTUs newVertexIndex leavesToMap choiceOpt
  | length leavesToAdd == nOTUs =
  let (eVertex, uVertex, edgeWeight) =  orderEdge $ getStartingPair choiceOpt distMatrix
      initialTree = (V.fromList[eVertex, uVertex],V.fromList [(eVertex, uVertex, edgeWeight)])
      leavesToAdd' = V.filter (/= eVertex) $ V.filter (/= uVertex) leavesToAdd
  in
  wagBest distMatrix initialTree leavesToAdd' nOTUs nOTUs leavesToAdd' choiceOpt
  | V.null leavesToAdd = (inTree, leavesToAdd, distMatrix)
  | otherwise =
  let addPosVect = V.map (addTaxonToTree distMatrix inTree leavesToAdd newVertexIndex) leavesToMap
    -- addPosVect = parmap rseq (addTaxonToTree distMatrix inTree leavesToAdd newVertexIndex) leavesToMap
      (firstLeafCost, _, _, _ ) = V.head addPosVect -- To initialize below
      (_, newTree, newLeavesToAdd, augmentedDistMatrix) = getBestLeafAdd (V.tail addPosVect) firstLeafCost (V.head addPosVect)
  in
  let progress = show  ((fromIntegral (100 * (newVertexIndex - nOTUs))/fromIntegral (nOTUs - 2)) :: Double)
  in
  trace (takeWhile (/='.') progress ++ "%")
  wagBest augmentedDistMatrix newTree newLeavesToAdd nOTUs (newVertexIndex + 1) newLeavesToAdd choiceOpt


-- | calculateWagnerTrees takes an input distance matrix (and options later) and returns
-- a tree (V,E) discription of Wagner tree with labelled internal veritices and branch lengths 
calculateWagnerTrees :: M.Matrix Double -> String -> (Tree, M.Matrix Double)
calculateWagnerTrees distMatrix choiceOpt =
  if M.dim distMatrix == (0,0) then error "Null distance matrix"
  else
    -- get initial pair of leaves and create initial tree
    let nOTUs = M.cols distMatrix
        allLeaves = V.fromList [0..(nOTUs - 1)]
    in
    let (newTree, _, endMatrix) = wagBest distMatrix (V.empty, V.empty) allLeaves nOTUs nOTUs allLeaves choiceOpt
    in
    (newTree, endMatrix)

-- | getEdges takes root Index and determines edges from root
getEdgeRoot :: Int -> V.Vector Edge -> Edge
getEdgeRoot edgeIndex edgeVect =
  if V.null edgeVect then (-1,-1,-1.0)
  else
   let (eVect, uVect, _) = V.head edgeVect
   in
   if (eVect == edgeIndex) || (uVect == edgeIndex) then V.head edgeVect else getEdgeRoot edgeIndex (V.tail edgeVect)

-- | getEdgesNonRoot takes root Index and determines edges from root returns String of Taxa
-- alwasy getting ordered edges and ordered tree
-- so eVect > uVect always; eVect can never be < nOTUs
--Need to add smaller tree left, bigger right
getEdgesNonRoot :: Int -> V.Vector Edge -> Int -> V.Vector String -> String
getEdgesNonRoot edgeIndex edgeVect nOTUs leafNames =
  --trace (show edgeIndex) (
  let terminal = (-1,-1,-1.0)
  in
  if V.null edgeVect then "END"
  else
   let thisEdge = getEdgeRoot edgeIndex edgeVect
   in
   if thisEdge == terminal then "ERROR"
   else
      let (eVect, uVect, weight) = thisEdge
          remainderEdges = V.filter (/= thisEdge) edgeVect
          eDesc = getEdgeRoot eVect remainderEdges
          uDesc = getEdgeRoot uVect remainderEdges
          eSubTree = getEdgesNonRoot eVect remainderEdges nOTUs leafNames
          uSubTree = getEdgesNonRoot uVect remainderEdges nOTUs leafNames
      in
      if V.null remainderEdges then
        (leafNames V.! uVect) ++ ":" ++ showDouble precision weight  ++ ","

      else if eVect == edgeIndex then

        if eVect < nOTUs then
          if uDesc /= terminal then "(" ++ (leafNames V.! eVect) ++ ":" ++ showDouble precision weight ++ "," ++ uSubTree ++ ")"
          else (leafNames V.! eVect) ++ ":" ++ showDouble precision weight ++ ","
        else
        if uVect < nOTUs then
          if eDesc /= terminal then "(" ++ (leafNames V.! uVect) ++ ":" ++ showDouble precision weight ++ "," ++ eSubTree ++  ")"
          else (leafNames V.! uVect) ++ ":" ++ showDouble precision weight ++ ","
        else
            if (eDesc /= terminal) && (uDesc == terminal) then eSubTree  ++ ":" ++ showDouble precision weight  ++ ","
            else if (eDesc == terminal) && (uDesc /= terminal) then uSubTree ++ ":" ++ showDouble precision weight ++ ","
            else
              if length eSubTree < length uSubTree then "(" ++ eSubTree ++ "," ++ uSubTree ++ ":" ++ showDouble precision weight ++ ")"
              else "(" ++ uSubTree ++ "," ++ eSubTree ++ ":" ++ showDouble precision weight ++ ")"

      else if uVect == edgeIndex then
        if uVect < nOTUs then
          if eDesc /= terminal then "(" ++ (leafNames V.! uVect) ++ ":" ++ showDouble precision weight ++  "," ++ eSubTree ++ ")"
          else (leafNames V.!uVect) ++ ":" ++ showDouble precision weight ++ ","

        else if eVect < nOTUs then
          if uDesc /= terminal then "(" ++ (leafNames V.! eVect) ++ ":" ++ showDouble precision weight ++ "," ++ uSubTree ++ ")"
          else (leafNames V.!eVect) ++ ":" ++ showDouble precision weight  ++ ","

        else
          if (eDesc /= terminal) && (uDesc == terminal) then eSubTree  ++ ":" ++ showDouble precision weight ++ ","
          else if (eDesc == terminal) && (uDesc /= terminal) then uSubTree ++ ":" ++ showDouble precision weight ++ ","
          else
              if length eSubTree < length uSubTree then "(" ++ eSubTree ++ "," ++ uSubTree ++ ":" ++ showDouble precision weight ++ ")"
              else "(" ++ uSubTree ++ ":" ++ showDouble precision weight ++ "," ++ eSubTree ++ ")"

      else getEdgesNonRoot edgeIndex remainderEdges nOTUs leafNames ++ ":" ++ showDouble precision weight ++ ","
    --))

-- | convertToNewick wrapper to remove double commas
convertToNewick :: V.Vector String ->Int ->  Tree -> String
convertToNewick leafNames outGroup wagTree = removeCrap $ convertToNewickGuts leafNames outGroup wagTree

-- | removeDoubleCommas removes second of double comas ",," -> "," 
-- this a hack to fix problem in convertToNewick
removeCrap :: String -> String
removeCrap inString =
  if length inString == 1 then inString
  else
    let firstChar = head inString
        secondChar = inString !! 1
    in
    if firstChar == ',' && secondChar == ',' then ',' : removeCrap (drop 2 inString)
    else if firstChar == ',' && secondChar == ')' then ')' : removeCrap (drop 2 inString)
    else firstChar : removeCrap (tail inString)


-- | convertToNewick converts Tree rep to Newick STring
-- includes edge cost--splits root edge cost into halves to help
-- tree viewers like FigTree
-- NEED TO ADD smaller group left larger group right for more legible trees
convertToNewickGuts :: V.Vector String ->Int ->  Tree -> String
convertToNewickGuts leafNames outGroup wagTree =
  let (inLeaves, inEdges) = wagTree
      newEdges = fmap orderEdge inEdges
      (_, edgeVect) = orderTree (inLeaves, newEdges)
      foundEdge = getEdgeRoot outGroup edgeVect
  in
  let (firstVert, secondVert, weight) = foundEdge
      remainderEdges = V.filter (/= foundEdge) edgeVect
  in
  -- this is embarassing bullshit  -- converting  ",,"  to ","
  if firstVert == outGroup then "(" ++ (leafNames V.! outGroup)  ++ ":" ++ showDouble 8 (weight/2.0) ++ "," ++ getEdgesNonRoot secondVert remainderEdges (V.length leafNames) leafNames ++ ":" ++ showDouble 8 (weight/2.0) ++ ")"
  else "(" ++ (leafNames V.! outGroup)  ++ ":" ++ showDouble 8 (weight/2.0) ++ "," ++ getEdgesNonRoot firstVert remainderEdges (V.length leafNames) leafNames ++ ":" ++ showDouble 8 (weight/2.0) ++ ")"
  --)
    --"();"

-- | getEdgeCost returns weight form edge tuple
getEdgeCost :: (Vertex, Vertex, Weight) -> Double
getEdgeCost (_, _, edgeWeight) = edgeWeight


-- | makeTreeFromOrder  takes an input order and other arguemnts and cretes tree using a single additoin
-- seqeunce, best plaecment for the leaf each round
makeTreeFromOrder :: M.Matrix Double -> Tree -> Int -> Int -> V.Vector Int -> (Tree, M.Matrix Double)
makeTreeFromOrder distMatrix initialTree nOTUs vertexIndex leavesToAdd =
  if null leavesToAdd then (initialTree, distMatrix)
  else
      let leaf = V.head leavesToAdd
          (_, newTree, _, augmentedDistMatrix) = addTaxonToTree distMatrix initialTree leavesToAdd vertexIndex leaf
      in
      makeTreeFromOrder augmentedDistMatrix newTree nOTUs (vertexIndex + 1) (V.tail leavesToAdd)

-- | getRandomAdditionSequence initializes based on input sequence and adds in order from there
getRandomAdditionSequence :: V.Vector String -> M.Matrix Double -> Int -> V.Vector Int -> TreeWithData
getRandomAdditionSequence leafNames distMatrix outgroup initiaLeavesToAdd =
  let nOTUs = V.length leafNames
  in
  let eVertex = initiaLeavesToAdd V.! 0
      uVertex = initiaLeavesToAdd V.! 1
      edgeWeight = distMatrix M.! (eVertex, uVertex)
      initialTree = (V.fromList[eVertex, uVertex],V.fromList [(eVertex, uVertex, edgeWeight)])
      leavesToAdd = V.filter (/= eVertex) $ V.filter (/= uVertex) initiaLeavesToAdd
  in
  let thisTree = makeTreeFromOrder distMatrix initialTree nOTUs nOTUs leavesToAdd
      -- (_, edgeVect) = fst thisTree
      treeCost = getTreeCost $ fst thisTree -- V.sum $ V.map getEdgeCost edgeVect
      newickTree = convertToNewick leafNames outgroup (fst thisTree) ++ "[" ++ showDouble precision treeCost ++ "]" ++ ";"
  in
  (newickTree, fst thisTree, treeCost, snd thisTree)

-- | doWagnerS takes user options and produces the Wagner tree methods desired (best, asis, or random)
-- outputs newick rep list 
doWagnerS :: V.Vector String -> M.Matrix Double -> String -> Int -> String -> [V.Vector Int]-> [TreeWithData]
doWagnerS leafNames distMatrix firstPairMethod outgroup addSequence replicateSequences =
  let nOTUs = V.length leafNames
  in
  if addSequence == "best" then
     let wagnerResult = calculateWagnerTrees distMatrix firstPairMethod
         -- (_, edgeVect) = fst wagnerResult
         treeCost = getTreeCost $ fst wagnerResult --- V.sum $ V.map getEdgeCost edgeVect
         newickTree = convertToNewick leafNames outgroup (fst wagnerResult) ++ "[" ++ showDouble precision treeCost ++ "]" ++ ";"
      in
      [(newickTree, fst wagnerResult, treeCost, snd wagnerResult)]
  else if addSequence == "asis" then
      let initialTree = (V.fromList[0, 1],V.fromList [(0, 1, distMatrix M.! (0,1))])
          leavesToAdd = V.fromList [2..(nOTUs-1)]
          asIsResult = makeTreeFromOrder distMatrix initialTree nOTUs nOTUs leavesToAdd
          treeCost = getTreeCost $ fst asIsResult -- V.sum $ V.map getEdgeCost asIsEdges
          newickTree = convertToNewick leafNames outgroup (fst asIsResult) ++ "[" ++ showDouble precision treeCost ++ "]" ++ ";"
      in
      [(newickTree, fst asIsResult, treeCost, snd asIsResult)]
  else if head addSequence == 'r' then
      if (length replicateSequences) == 0 then error "Zero replicate additions specified--could be error in configuration file"
      else 
        let (chunkSize, _) = quotRem (length replicateSequences) getNumThreads
            randomAddTrees = fmap (getRandomAdditionSequence leafNames distMatrix outgroup) replicateSequences `using` parListChunk chunkSize myStrategy -- was rseq not sure whats better
            -- randomAddTrees = parmap rseq (getRandomAdditionSequence leafNames distMatrix outgroup) replicateSequences
        in
        randomAddTrees
  else error ("Addition sequence " ++ addSequence ++ " not implemented")

-- | getRandomReps takes teh buyild arguemnrt and returns the number of random replicates 0 if not ranomd the numbewr otherwise
getRandomReps :: String -> Int
getRandomReps inString
  | head inString /= 'r' = 0
  | length inString < 8 = error ("Need to specify repicate number after '=' in " ++ inString)
  | otherwise = read (drop 7 inString) :: Int


-- | ranList generates random list of positive integers
ranList :: Rand.StdGen -> Int -> Int -> [Int]
ranList sg n maxValue = take n $ Rand.randomRs (0,maxValue) sg

-- | driver function to generate list of positive integers 
{-# NOINLINE betterRandomList #-}
betterRandomList :: Int -> Int -> [Int]
betterRandomList n maxValue= unsafePerformIO $ do
    sg <- Rand.getStdGen
    return $ ranList sg n maxValue

-- | getNumThreads gets number of COncurrent  threads
{-# NOINLINE getNumThreads #-}
getNumThreads :: Int
getNumThreads = unsafePerformIO getNumCapabilities


-- | getRandAddVect takes list of sequence integers, randomizes and returns as vector
{-# NOINLINE getRandAddVect #-}
getRandAddVect :: [a] -> V.Vector a
getRandAddVect inList =
  let randIntList = betterRandomList (length inList) (length inList -1)
  in
  trace (show randIntList)
  V.fromList $ RandS.shuffle inList randIntList


-- getBestTrees takes newick and sorts on comment at end with cost
getBestTrees :: String -> Int -> [TreeWithData] -> Double -> [TreeWithData] -> [TreeWithData]
getBestTrees keepMethod number inList curBestCost curBestTrees =
  if null inList then
      -- apply keep method if not keeping all
      if length curBestTrees <= number then curBestTrees
      else if keepMethod == "first" then take number curBestTrees
      else if keepMethod == "last" then reverse $ take number $ reverse curBestTrees
      else if keepMethod == "random" then
        -- not Fitch but shuffles elements and takes first n
        let randIntList = betterRandomList (length curBestTrees) (length curBestTrees - 1)
            newList =  RandS.shuffle curBestTrees randIntList
        in
        take number newList
      else error ("Keep method " ++ keepMethod ++ " not implemented")
  else
    let firstTree = head inList
        (_, _, firstCost, _) = firstTree
    in
    if firstCost < curBestCost then getBestTrees keepMethod number (tail inList) firstCost [firstTree]
    else if withinEpsilon firstCost curBestCost then getBestTrees keepMethod number (tail inList) firstCost (firstTree : curBestTrees)
    else getBestTrees keepMethod number (tail inList) curBestCost curBestTrees


-- | getUniqueTrees saves uniqe newick trees
-- different paths--ehnce different distMatrices coud result in same newick tree
getUniqueTrees :: [TreeWithData] -> [TreeWithData] -> [TreeWithData]
getUniqueTrees inList uniqueList =
  if null inList then uniqueList
  else
    let firstTree = head inList
        fstNewick = fst4 firstTree
        -- (firstTreeNewick, firstTreeTree) = firstTree
    in
    if fstNewick `notElem` (fmap fst4 uniqueList) then getUniqueTrees (tail inList) (firstTree : uniqueList)
    else getUniqueTrees (tail inList) uniqueList

-- | keepTrees filters newick trees based on options
-- all keep all
-- best shortest (and unique) allows number of max to save 
-- unique unique  representations irespective of length
-- keep metyhod for save first | last | atRandom if buffer full
keepTrees :: [TreeWithData] -> String -> String -> Double -> [TreeWithData]
keepTrees inList saveMethod keepMethod curBestCost
  | null inList = []
  | saveMethod == "all" = inList
  | take 6 saveMethod == "unique" =
    if length saveMethod == 6 then getUniqueTrees inList []
    else
      let number = read (drop 7 saveMethod) :: Int
      in
      take number $ sortOn thd4 $ getUniqueTrees inList []
  | take 4 saveMethod == "best" =
    if length saveMethod == 4 then getUniqueTrees (getBestTrees keepMethod (maxBound :: Int) inList curBestCost []) []
    else if (saveMethod !! 4) == ':' then
      let number =  read (drop 5 saveMethod) :: Int
          saveTrees = take number $ sortOn thd4 $ getUniqueTrees inList []
          (_, _, bestCost, _) = head saveTrees
      in
      getBestTrees keepMethod number saveTrees bestCost []
    else error ("Save method " ++ saveMethod ++ " improperly formatted")
  | otherwise = error ("Save method " ++ saveMethod ++ " not implemented")

-- | edgeHasVertex takes an vertex and an edge and returns Maybe Int
-- of other vertex
edgeHasVertex :: Vertex -> Edge -> Maybe (Vertex, Edge)
edgeHasVertex inVert inEdge =
  let (a, b, _) = inEdge
  in
  if a == inVert then Just (b, inEdge)
  else if b == inVert then Just (a, inEdge)
  else Nothing

-- | getSubEdges take a Vector of edges and retuns a list of edges connected to input vertex
-- uses nOTUs to know when to stop recursing
getSubEdges :: [Vertex] -> Int -> V.Vector Edge -> V.Vector Edge -> String -> V.Vector Edge
getSubEdges inVertexList nOTUs edgeVect subEdgeVect howMany
  | V.null edgeVect = subEdgeVect
  | null inVertexList = subEdgeVect
  | otherwise =
    let startVertex = head inVertexList
        foundVect = V.filter (/= Nothing) $ V.map (edgeHasVertex startVertex) edgeVect
    in
    if V.null foundVect then getSubEdges (tail inVertexList) nOTUs edgeVect subEdgeVect howMany-- only terminals left
    else if V.length foundVect /= 2 then error ("Index (" ++ howMany ++ ")" ++ show startVertex ++ "->found " ++ show (V.length foundVect) ++ " but should be two edges in " ++ show foundVect ++ " in " ++ show edgeVect)
    else
      let thingsFound = V.map fromJust foundVect
          verticesFound = V.toList $ V.map fst thingsFound
          edgesFound = V.map snd thingsFound

          -- add found edges to subEdgeSet
          newSubEdgeVect =  subEdgeVect V.++ edgesFound
          -- delete delete edges from those to be searched
          newEdgeVect = V.filter (/= V.last edgesFound) $ V.filter (/= V.head edgesFound) edgeVect
      in
      -- recurse on vertices that were found
      if howMany == "first2" then edgesFound
      else getSubEdges (verticesFound ++ tail inVertexList) nOTUs newEdgeVect newSubEdgeVect howMany

-- | Subtrace vector subtracts elements of vector a from vector b
-- is thins n^2 ?
-- edges are directed
subtractVector :: V.Vector Edge -> V.Vector Edge -> V.Vector Edge
subtractVector a b
  | V.null a = b
  | V.null b = V.empty
  | otherwise =
    let firstB = V.head b
        notFound = V.notElem firstB a
    in
    if notFound then V.cons firstB (subtractVector a (V.tail b))
    else subtractVector a (V.tail b)


-- | adjustInternalEdgeVertex adjusts vertex of an internal (ie non-pendent, none 2 terminal edge)
-- assumes hv > lv
adjustInternalEdgeVertex :: Vertex -> Vertex -> Vertex -> Int -> Int -> Vertex
adjustInternalEdgeVertex inV hV lV maxOffSet nOTUs
  | inV <= max nOTUs lV = inV
  | inV > hV = inV - maxOffSet
  | inV > lV = inV - 1
  | otherwise = error ("This can't happen " ++ show (inV, lV, hV))
  -- )

-- | adjustVertex reduces vertex (hV,lV,_) is the dge that was deleted
-- index by 2,1, or 0 if greater than hVert, lVert, or neither
-- asumes hv > lv; and inE > inU
-- if selfe edge  its a termin and only reduce by 1
-- assumes edge is ordered e > u
adjustVertex :: Edge -> Vertex -> Vertex -> Int -> Edge
adjustVertex (inE, inU, w) hV lV nOTUs
  | (inE <= max lV nOTUs) && (inU <= max lV nOTUs) = (inE, inU, w)
  | lV < nOTUs =  -- update pendant edge was deleted since hV > lV; hV must be > nOTUs   
    if (inE > max hV nOTUs) && (inU > max hV nOTUs) then (inE - 1, inU - 1, w)
    else if (inE <= max hV nOTUs) && (inU > max hV nOTUs) then (inE, inU - 1, w)
    else if (inE > max hV nOTUs) && (inU <= max hV nOTUs) then (inE - 1, inU, w)
    else (inE, inU, w)
  | otherwise =                                                          -- internal edge was deleted                                                         -- both deleted verteces internal
    let newE = adjustInternalEdgeVertex inE hV lV 2 nOTUs
        newU = adjustInternalEdgeVertex inU hV lV 2 nOTUs
    in
    (newE, newU, w)

-- | updateVertexNUmbersOnEdges taked vertex numbers and updates HTU indices to reflect 
-- the deletion of those vertices 
-- ASSUMES edges are ordered (a,b,weight) a > b 
-- subtracts 2 from index if > the bigger of teh two, subtract 1 if bigger than lower,
-- otherwise leaves unchanged
updateVertexNUmbersOnEdges  :: Vertex -> Vertex -> V.Vector Edge -> Int -> V.Vector Edge
updateVertexNUmbersOnEdges eVert uVert edgeList nOTUs =
  -- trace ("In updateVertexNUmbersOnEdges") (
  if V.null edgeList then V.empty
  else
      let hVert = max eVert uVert
          lVert = min eVert uVert
          newVertex = adjustVertex (V.head edgeList) hVert lVert nOTUs
      in
      V.cons newVertex (updateVertexNUmbersOnEdges eVert uVert (V.tail edgeList) nOTUs)
     -- )

{-
-- | updateRow updates row of distance matrix (shortens) based on deleteing two vertex entries
-- assume 1st is higher and second lower
updateRow :: Vertex -> Vertex -> [Double] -> Int -> Int -> Int -> Edge -> Edge -> [Double]
updateRow hVert lVert rowList nOTUs rowCounter columnCounter c1Edge c2Edge =
  let (c1E, c1U, c1W) = c1Edge
      (c2E, c2U, c2W) = c2Edge
  in
  if null rowList then []
  else
      if columnCounter < nOTUs then
          if ((rowCounter == c1E) && (columnCounter == c1U)) || ((rowCounter == c1U) && (columnCounter == c1E)) then c1W : updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1) c1Edge c2Edge else (if ((rowCounter == c2E) && (columnCounter == c2U)) || ((rowCounter == c2U) && (columnCounter == c2E)) then c2W : updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1) c1Edge c2Edge else head rowList : updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1) c1Edge c2Edge)
      else if (columnCounter == hVert) || (columnCounter == lVert) then updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1) c1Edge c2Edge
      else
          if ((rowCounter == c1E) && (columnCounter == c1U)) || ((rowCounter == c1U) && (columnCounter == c1E)) then c1W : updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1) c1Edge c2Edge else (if ((rowCounter == c2E) && (columnCounter == c2U)) || ((rowCounter == c2U) && (columnCounter == c2E)) then c2W : updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1) c1Edge c2Edge else head rowList : updateRow hVert lVert (tail rowList) nOTUs rowCounter (columnCounter + 1)  c1Edge c2Edge)


-- | makeNewRows updates rows recursively based on deleting rows and columns in hVert (greater)  and lVert (lesser)
makeNewRows :: Vertex -> Vertex -> [[Double]]  -> Int -> Int -> Edge -> Edge -> [[Double]]
makeNewRows hVert lVert distListList nOTUs rowCounter c1Edge c2Edge
  | null distListList = []
  | rowCounter < nOTUs =
  let firstRow = head distListList
      newRow = updateRow hVert lVert firstRow nOTUs rowCounter 0 c1Edge c2Edge
  in
  newRow : makeNewRows hVert lVert (tail distListList) nOTUs (rowCounter + 1) c1Edge c2Edge
  | (rowCounter == hVert) || (rowCounter == lVert) = makeNewRows hVert lVert (tail distListList) nOTUs (rowCounter + 1) c1Edge c2Edge
  | otherwise =
  let firstRow = head distListList
      newRow = updateRow hVert lVert firstRow nOTUs rowCounter 0 c1Edge c2Edge
  in
  newRow : makeNewRows hVert lVert (tail distListList) nOTUs (rowCounter + 1) c1Edge c2Edge


-- | updateDistMatrix updates distance matrix to remove eVertex and uVertex columns and rows as in updateVertexNUmbersOnEdges
--  update costs for two contracte edges c1Edge and c2Edge
updateDistMatrix :: Vertex -> Vertex ->  M.Matrix Double -> Int -> Edge -> Edge -> M.Matrix Double
updateDistMatrix eVert uVert distMatrix nOTUs c1Edge c2Edge =
  let hVert = max eVert uVert
      lVert = min eVert uVert
      distListList = M.toLists distMatrix
  in
  M.fromLists $ makeNewRows hVert lVert distListList nOTUs 0 c1Edge c2Edge
-}

-- | updateDistMatrix updates distance matrix to remove eVertex and uVertex columns and rows as in updateVertexNUmbersOnEdges
-- update costs for two contracte edges c1Edge and c2Edge
-- the distance update takes place first and row/coumn removal second
-- this to keep the proper indices (since the non-deleted rows and columns indices are changed) 
updateDistMatrix :: Vertex -> Vertex ->  M.Matrix Double -> Int -> Edge -> Edge -> M.Matrix Double
updateDistMatrix eVert uVert distMatrix nOTUs c1Edge c2Edge =
  let validEdgeList = filter ((>= 0).fst3) [c1Edge, c2Edge]
      newMatrix = M.unsafeUpdateMatrix distMatrix validEdgeList
      newMatrix' = M.deleteRowsAndColumns newMatrix (filter (> (nOTUs - 1)) [eVert, uVert])
  in
  -- trace (M.showMatrixNicely newMatrix') 
  newMatrix'
  


-- | getEndVertices takes a pair of edges and returns the non-index vertices
getEndVertices :: V.Vector Edge -> Vertex -> (Vertex, Vertex)
getEndVertices inEdges index =
  if V.length inEdges /= 2 then error ("Edge number should be 2 not " ++ show (V.length inEdges) ++ " in getEndVertices")
  else
    let (fEVertex, fUVertex, _) = V.head inEdges
        (sEVertex, sUVertex, _) = V.last inEdges
    in
    if fEVertex == index then
      if sEVertex == index then (fUVertex, sUVertex)
      else
          (fUVertex, sEVertex)
    else if fUVertex == index then
      if sEVertex == index then
         (fEVertex, sUVertex)
      else
          (fEVertex, sEVertex)
    else error ("Error findeing ends of " ++ show (V.head inEdges) ++ " and " ++ show (V.last inEdges))

-- | getContractedWeight gets edge contractd edge weights ither form distMatrix for  OTUs or 1 OTU and 1 
-- internal wertex or by 4 point metric (maximum of two estimations)
getContractedWeight :: Vertex -> Vertex -> M.Matrix Double -> Int -> V.Vector Edge -> Double
getContractedWeight aVert bVert distMatrix nOTUs edgeVect
  | aVert < nOTUs && bVert < nOTUs = distMatrix M.! (aVert, bVert)
  | aVert < nOTUs = distMatrix M.! (aVert, bVert)
  | bVert < nOTUs = distMatrix M.! (aVert, bVert)
  | otherwise =
  -- both internal 4-point mertic estimation
  --get edges connected to the contracted edge
  let aEdgeVect = getSubEdges [aVert] nOTUs edgeVect V.empty "first2"
      bEdgeVect = getSubEdges [bVert] nOTUs edgeVect V.empty "first2"
      (a,b) = getEndVertices aEdgeVect aVert
      (c,d) = getEndVertices bEdgeVect bVert
      firstEstimate  = (distMatrix M.! (a,c)) + (distMatrix M.! (b,d)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
      secondEstimate = (distMatrix M.! (a,b)) + (distMatrix M.! (b,c)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
  in
  max firstEstimate secondEstimate / 2.0


-- | contractEdges takes two edges that share a vertex and fuses them
-- If terminal then return selfdge
-- update contracted edge weight
-- dummy edge (-1,-1,-1) is returned when no edge to contract so won't be updated in distance Matrix
contractEdges :: M.Matrix Double -> Int -> V.Vector Edge -> Vertex -> V.Vector Edge -> Edge
contractEdges distMatrix nOTUs edgeVect index allEdges
  | V.null edgeVect = (-1,-1,0)
  | V.length edgeVect == 1 =
    let (a,b,w) = V.head edgeVect
    in
    if (a == index) || (b == index) then (index, index, w)
    else error ("Contacting single edge: " ++ show (V.head edgeVect) ++ " with vertex " ++ show index ++ " not found")
  | otherwise =
    let (aVert,bVert) = getEndVertices edgeVect index
        newWeight = getContractedWeight aVert bVert distMatrix nOTUs (subtractVector edgeVect allEdges)
    in
    orderEdge (aVert, bVert, newWeight)

-- | getDistance creates distances from remaing OTUs to new vertex (rowID) via Farris 1972
-- rowID is the row (or column) of distance Matrix
getDistance :: M.Matrix Double -> Double -> Double -> Double -> Int -> Int -> Int -> Int -> Double
getDistance origDist addCost eVertLeafDist uVertLeafDist leafIndex eVertex uVertex rowID =
  let first  = (origDist M.! (rowID, leafIndex)) - addCost
      second = (origDist M.! (rowID, eVertex)) - eVertLeafDist
      third  = (origDist M.! (rowID, uVertex)) - uVertLeafDist
  in
  maximum [first, second, third]

-- | getNewDistMatrix takes distMatrix and adds Cost for new eVertLeafDist uVertLeafDist
-- created in build process (new HTUs)
-- should be complete (on input) for leaves already added (initail paiwise distances and HTUs added) 
-- adds a single new row (minus last 0.0 as new row at end) which is appended
getNewDistMatrix :: M.Matrix Double -> Double -> Double -> Double -> Int -> Int -> Int -> M.Matrix Double
getNewDistMatrix origDist addCost eVertLeafDist uVertLeafDist eVertex uVertex leafIndex =
    let columnHolder = V.fromList [0..(M.rows origDist - 1)] -- List of HTU and OTU indices in pairwise dist matrix
        newDistRow = V.map (getDistance origDist addCost eVertLeafDist uVertLeafDist leafIndex eVertex uVertex) columnHolder
        newDistRow' = newDistRow `V.snoc` (0.0 :: Double) 
    in
    M.addMatrixRow origDist newDistRow'
    {-
    --This is a waste--change structure
    let oldRows = V.fromList (M.toRows origDist)
        newRows = addNewColumn oldRows newDistRow newDistRow
    in
    M.fromRows (V.toList newRows)
    --)
    -}

-- | enterNewEdgeCost cretes new row of costs from edge vectors 
-- infty NT.infinity if not there
-- this makes update n^2 which is dumb
enterNewEdgeCost :: Int -> V.Vector Edge -> Int -> Double
enterNewEdgeCost columnNumber edgeVect rowNumber =
  if V.null edgeVect then NT.infinity --not found so Infty
  else
    let (a, b, weight) = V.head edgeVect
    in
    if (columnNumber == a && rowNumber == b) || (columnNumber == b && rowNumber == a) then weight else enterNewEdgeCost columnNumber (V.tail edgeVect) rowNumber

-- | addEdgesToDistMatrix adds new edges from internal node join to existing matrix
-- 0 if not created so only the two new vertices and the edges they touch (3 each)
-- so need to add two columns and two rows, for new vertices I, and II (here numIn and numIn + 1)
getNewDistMatrixInternal :: M.Matrix Double -> V.Vector Edge -> M.Matrix Double
getNewDistMatrixInternal inMatrix newEdgeVect =
  -- trace ("In getNewDistMatrixInternal") (
  if V.length newEdgeVect /= 5 then error ("Wrong size edgeVector shoud be 5 and is " ++ show (V.length newEdgeVect))
  else
    let numIn = M.rows inMatrix
        columnHolder = [0..(numIn - 1)]
        -- newDistColumnI = fmap (enterNewEdgeCost numIn newEdgeVect) columnHolder ++ [0.0, enterNewEdgeCost numIn newEdgeVect (numIn + 1)]
        -- newDistColumnII = fmap (enterNewEdgeCost (numIn + 1) newEdgeVect) columnHolder ++ [enterNewEdgeCost numIn newEdgeVect (numIn + 1), 0.0]
        newDistColumnI = (fmap (enterNewEdgeCost numIn newEdgeVect) columnHolder) ++ [0.0]
        newDistColumnII = (fmap (enterNewEdgeCost (numIn + 1) newEdgeVect) columnHolder) ++ [enterNewEdgeCost numIn newEdgeVect (numIn + 1), 0.0]
        {-
        --This is a waste--change structure
        oldRows = M.toLists inMatrix
        newRows = addTwoColumns oldRows newDistColumnII newDistColumnII
        newLists = newRows ++ [newDistColumnI, newDistColumnII]
        -}
    in
    M.addMatrices inMatrix  (V.fromList [V.fromList newDistColumnI, V.fromList newDistColumnII])
    {-
    M.fromLists newLists
    -- )
    -}

-- | connectEdges takes two vectors of edges and adds and edge between the two in edgesToConnect
-- this deletes the two old edges from the edge Vectors and creates a new tree with the five new edges
-- added in; the addition cost is for total of created edges minus edges that were destroyed
connectEdges :: M.Matrix Double -> V.Vector Edge -> V.Vector Edge -> V.Vector Edge -> (Double, Tree, M.Matrix Double)
connectEdges distMatrix eEdges uEdges edgesToConnect
  | V.null eEdges = error "Empty e Edge vector in connectEdges"
  | V.null uEdges = error "Empty u Edge vector in connectEdges"
  | V.length edgesToConnect /= 2 = error ("There shoule be 2 edges to connect and ther are " ++ show (V.length edgesToConnect))
  | otherwise =
  let edgesKept = subtractVector edgesToConnect eEdges V.++ subtractVector edgesToConnect uEdges
      numInTree = M.rows distMatrix
      -- order new edges 
      (a, b, wAB) = V.head edgesToConnect
      (c, d, wCD) = V.last edgesToConnect
      firstEstimate  = (distMatrix M.! (a,c)) + (distMatrix M.! (b,d)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
      secondEstimate = (distMatrix M.! (a,b)) + (distMatrix M.! (b,c)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
      centralEdgeCost  = max firstEstimate secondEstimate / 2.0
      centralEdge = orderEdge (numInTree, numInTree + 1, centralEdgeCost)
      aCost = ((distMatrix M.! (c,a)) - (distMatrix M.! (c,b)) + (distMatrix M.! (a,b))) / 2.0
      bCost = (distMatrix M.! (a,b)) - aCost
      cCost = ((distMatrix M.! (a,c)) - (distMatrix M.! (a,d)) + (distMatrix M.! (c,d))) / 2.0
      dCost = (distMatrix M.! (c,d)) - cCost
      newAEdge = orderEdge (a, numInTree, aCost) -- edge cost not unique either (a or b, numInTree) could be max or min
      newBEdge = orderEdge (b, numInTree, bCost) -- edge cost not unique either (a or b, numInTree) could be max or min
      newCEdge = orderEdge (c, numInTree + 1, cCost) -- edge cost not unique either (a or b, numInTree + 1) could be max or min
      newDEdge = orderEdge (d, numInTree + 1, dCost) -- edge cost not unique either (a or b, numInTree + 1) could be max or min
      newEdgeVect = V.fromList [centralEdge, newAEdge, newBEdge, newCEdge, newDEdge]
      newDistMatrix = getNewDistMatrixInternal distMatrix newEdgeVect
  in
  (centralEdgeCost + aCost + bCost + cCost + dCost - wAB - wCD, (V.empty, newEdgeVect V.++ edgesKept), newDistMatrix)


-- | addEdgeToSplit adds a new edge to a specifc pair of input edges detemined addition cost
-- and creted new edge set with appropriate weights
addEdgeToSplit :: V.Vector Edge -> V.Vector Edge -> Edge -> Edge -> M.Matrix Double -> V.Vector Edge -> (Double, Tree, M.Matrix Double)
addEdgeToSplit eEdges uEdges eTerminal uTerminal distMatrix edgesToConnect
  | V.null eEdges &&  V.null uEdges = error "Empty e/u Edge vectors in addEdgeToSplit"
  | V.null edgesToConnect = error "Empty eEdge vector in addEdgeToSplit"
  | fst3 (V.head eEdges) == (-1) = -- adding eOTU, V.last since the (-1,-1,0) edge should always be first
      addToEdgeSwap distMatrix (fst3 eTerminal) (V.empty, uEdges) (M.rows distMatrix) (V.last edgesToConnect)
  | fst3 (V.head uEdges) == (-1) = -- adding uOTU
      addToEdgeSwap distMatrix (fst3 uTerminal) (V.empty, eEdges) (M.rows distMatrix) (V.last edgesToConnect)
  | otherwise = -- both internal edges
      connectEdges distMatrix eEdges uEdges edgesToConnect


-- | getVertexSet take a vector of edges and creates the set of vertex numbers
getVertexSet :: V.Vector Edge -> Set.Set Vertex
getVertexSet edgeVect =
  if V.null edgeVect then Set.empty
  else
      let (a,b,_) = V.head edgeVect
          thisSet =  Set.fromList [a,b]
      in
      Set.union thisSet (getVertexSet $ V.tail edgeVect)

-- | edgeHasVertexBool uses edgeHasVertex but returns True if found, False if not
edgeHasVertexBool :: Vertex -> Edge -> Bool
edgeHasVertexBool inVertex inEdge =
  isJust (edgeHasVertex inVertex inEdge)

-- | checkDegreeVerts checks degree if each vertex by number of edges containing that number in edgeVect
-- if < nOTUS then should be 1 otherwise 3
checkDegreeVerts :: Int -> Int -> V.Vector Edge -> String -> Int -> Bool
checkDegreeVerts nOTUs curVertex edgeVect whereString adjustment  =
  if curVertex == (2 * nOTUs) - 2 - adjustment then trace whereString True
  else
    let vertDegree = V.length $ V.findIndices (edgeHasVertexBool curVertex) edgeVect
    in
    if curVertex < nOTUs then -- leaf
      if vertDegree /= 1 then error (whereString ++ " Vertex " ++ show curVertex ++ " has degreee " ++ show vertDegree ++ " and should have 1 in edge set " ++ show edgeVect)
      else checkDegreeVerts nOTUs (curVertex + 1) edgeVect whereString adjustment
    else
      if vertDegree /= 3 then error (whereString ++ "Vertex " ++ show curVertex ++ " has degreee " ++ show vertDegree ++ " and should have 3 in edge set " ++ show edgeVect)
      else checkDegreeVerts nOTUs (curVertex + 1) edgeVect whereString adjustment

-- | isProperTree takes a set of edges and checks for connectedness (3 for each non-OTU vertex edge, 1 for each vertex) 
-- not yet checking for reaching all OTUs and HTUs (connected)
-- Does NOT check for vertex set; only an edge check
isProperTree :: Tree -> Int -> Bool
isProperTree inTree matDim =
  let (_, edgeVect) = inTree
      totalVerts = matDim
      nOTUs = (matDim + 2) `div` 2
      totalEdges = (2 * nOTUs) - 3
      vertexSet = getVertexSet edgeVect
  in
  if totalEdges /= V.length edgeVect then error ("Matrix implies " ++ show totalEdges ++ " edges " ++ " but there are " ++ show (V.length edgeVect))
  else if Set.size vertexSet /= totalVerts then error ("Matrix implies " ++ show totalVerts ++ " vertices " ++ " but there are " ++ show (Set.size vertexSet))
  else
      let terminalVertSet = Set.filter (< nOTUs) vertexSet
      in
      if Set.size terminalVertSet /= nOTUs then error ("Matrix implies " ++ show nOTUs ++ " OTUs " ++ " but there are " ++ show (Set.size terminalVertSet) ++ "vertex indices < nOTUs")
      else checkDegreeVerts nOTUs 0 edgeVect "Full check" 0 || error "Error in degree of vertices"

-- | debug dummy edge
getDummyVertex :: Int -> Int -> Int -> Int
getDummyVertex eV uV nOTUs =
    let first = min  eV uV
    in
    if first < nOTUs then first
    else (-1)

-- | splitTree takes a tree description and its edgeList and return pairs of edge list
-- split at input edge in tree with "repaird"/contracted edges, delta, and original
-- edge (pairs of vertices and weight) for each split
splitTree :: M.Matrix Double -> Tree -> Double -> Edge -> SplitTreeData
splitTree distMatrix inTree inTreeCost edgeToRemove =
  -- check if proper tree--remove later
  let (_, edgeVect) = inTree
      (eVertex, uVertex, _) = edgeToRemove

      -- newEdgeSet = subtractVector (V.cons edgeToRemove $ eEdges V.++ uEdges) edgeVect   
      newEdgeSet = V.filter (/= edgeToRemove) edgeVect
      nOTUs = div (3 + V.length edgeVect) 2
      eSubEdges = getSubEdges [eVertex] nOTUs newEdgeSet V.empty "all"
      uSubEdges = getSubEdges [uVertex] nOTUs newEdgeSet V.empty "all"

      -- get edges that need to be contracted and re-estimate weights 
      eEdges = getSubEdges [eVertex] nOTUs eSubEdges V.empty "first2"
      uEdges = getSubEdges [uVertex] nOTUs uSubEdges V.empty "first2"
      eMergedEdge = contractEdges distMatrix nOTUs eEdges eVertex edgeVect
      uMergedEdge = contractEdges distMatrix nOTUs uEdges uVertex edgeVect

      -- remove non-contracted edges and add in contracted edges
      eSubEdges' = V.cons eMergedEdge (subtractVector eEdges eSubEdges)
      uSubEdges' = V.cons uMergedEdge (subtractVector uEdges uSubEdges)

      -- need to know this order for SPR/TBR so which edge was in which set
      previousEdges = V.fromList [eMergedEdge, uMergedEdge]

      -- map new HTU indices in edges and create new distance matrix
      -- HTU indices are updated first--then values updated in distMatrix as rows/columns are
      -- deleted to remove info from delted edge
      eSubEdges'' = updateVertexNUmbersOnEdges eVertex uVertex (V.map orderEdge eSubEdges') nOTUs
      uSubEdges'' = updateVertexNUmbersOnEdges eVertex uVertex (V.map orderEdge uSubEdges') nOTUs
      previousEdges'' = updateVertexNUmbersOnEdges eVertex uVertex (V.map orderEdge previousEdges) nOTUs

      -- Update with deleted node and reestimated contracted edges
      -- update matrix the costs (2 ij x 2 ji) of contracted edges 
      distMatrix'' = updateDistMatrix eVertex uVertex distMatrix nOTUs eMergedEdge uMergedEdge -- (V.head previousEdges'') (V.last previousEdges'')

      -- Delta calcualted by differene in original tree cost and split and readdition to split edges (then tail splits later
      -- so not remake the original tree)
      splitCost = getTreeCost (V.empty, eSubEdges'') + getTreeCost (V.empty, uSubEdges'')

      -- Delta of tree length will be weights of the the edges removed - the weights of the two edges contracted and reestimated
      -- delta = weight + (V.sum $ V.map thd3 eEdges) +  (V.sum $ V.map thd3 uEdges) - (V.sum $ V.map thd3 previousEdges'')
      delta =  inTreeCost - splitCost   -- readditionCost
  in
  if eVertex < nOTUs then (V.singleton (eVertex, eVertex, 0.0), uSubEdges'', delta, previousEdges'', distMatrix'') --this so know a pendant edge
  else if uVertex < nOTUs then (V.singleton (uVertex, uVertex, 0.0), eSubEdges'', delta, previousEdges'', distMatrix'') --this so know a pendant edge
  else
      -- neded to have e first then u for SPR/TBR so can coordinate with previous edges
      (eSubEdges'', uSubEdges'', delta, previousEdges'', distMatrix'')


-- | sieveTrees takes a list of (addition cost, Tree, distnce matrix) and returns list
-- of better or equal trees to input delta
sieveTrees :: Double -> Double -> V.Vector (Double, Tree, M.Matrix Double) -> V.Vector String -> Int -> [TreeWithData] -> [TreeWithData]
sieveTrees inDelta curBestCost inAddList leafNames outgroup savedTrees =
  if null inAddList then savedTrees
  else
      let firstTuple = V.head inAddList
          (firstDelta, firstTree, firstMatrix) = firstTuple
          newCost = curBestCost - inDelta + firstDelta
          -- checkCost = getTreeCost firstTree
          newickTree = convertToNewick leafNames outgroup firstTree ++ "[" ++ show newCost ++ "]" ++ ";"
          newTuple = (newickTree, firstTree, newCost, firstMatrix)
      in
      if firstDelta > inDelta then sieveTrees inDelta curBestCost  (V.tail inAddList) leafNames outgroup savedTrees
      else
        if newCost < curBestCost then
            sieveTrees inDelta curBestCost (V.tail inAddList) leafNames outgroup [newTuple]
        else if withinEpsilon newCost curBestCost then sieveTrees inDelta curBestCost (V.tail inAddList) leafNames outgroup (newTuple : savedTrees)
        else sieveTrees inDelta curBestCost (V.tail inAddList) leafNames outgroup savedTrees


-- | reAddTerminals checks to see if split on pendant edge--if so reads terminal to each edge but saves equal cost if found
-- and saveMethod specifies it.  Idenitical to the wagner addition process with saveing equal as option 
-- could use addEdgeToSplit for consistancey with SPR/TBR
reAddTerminals :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> [TreeWithData]
reAddTerminals rejoinType curBestCost leafNames outGroup split =
  if rejoinType /= "otu" then error ("Incorrect swap function in reAddTerminals: " ++ rejoinType)
  else
    let (eEdgeVect, uEdgeVect, delta, _, distMatrix) = split
        nOTUs = V.length leafNames
    in
    if (V.length eEdgeVect > 1) || ((fst3 (V.head eEdgeVect) /= snd3 (V.head eEdgeVect)) && (fst3 (V.head eEdgeVect) < nOTUs) && (snd3 (V.head eEdgeVect) < nOTUs)) then []
    else (if M.rows distMatrix /= ((2 * nOTUs) - 3) then error ("Dist Matrix incorrect size " ++ show (M.dim distMatrix) ++ " should be " ++ show ((2 * nOTUs) - 3, (2 * nOTUs) - 3))
    else
      let newLeafIndex = M.rows distMatrix
      -- take tail of uEdgeVect so not regerate input tree 
          additionList = V.map (addToEdgeSwap distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex) uEdgeVect -- (V.tail uEdgeVect) -- tail so not hit original tree, leave all to reestimate if necesary
          minAdditionCost = V.minimum (V.map fst3 additionList)
      in
      if minAdditionCost > delta then []
      else sieveTrees delta curBestCost additionList leafNames outGroup [])


-- | add leaf to an edge creating new tree with distances and add cost, also augmented distance matrix
-- but this for swap so returns entire new3-edge cost  so not Farris triangle it is sum of three diveded by 2
addToEdgeSwapRecurse :: Double -> M.Matrix Double -> Int -> Tree -> Int -> V.Vector Edge -> (Double, Tree, M.Matrix Double)
addToEdgeSwapRecurse inDelta distMatrix leaf initialTree newLeafIndex inEdgeVect =
  if V.null inEdgeVect then (inDelta, initialTree, distMatrix)
  else 
    let inEdge@(eVertex, uVertex, inWeight) = V.head inEdgeVect
        (initialVertexVect, initialEdgeVect) = initialTree
        addCost = ((distMatrix M.! (leaf, eVertex)) + (distMatrix M.! (leaf, uVertex)) - (distMatrix M.! (eVertex, uVertex))) / 2.0
        eVertLeafDist = (distMatrix M.! (leaf, eVertex)) - addCost
        uVertLeafDist = (distMatrix M.! (leaf, uVertex)) - addCost
        newVertexVect = V.snoc initialVertexVect leaf
        newEdges = V.fromList [(leaf,newLeafIndex, addCost),(eVertex, newLeafIndex, eVertLeafDist),(uVertex, newLeafIndex, uVertLeafDist)]
        cleanupEdges = V.filter (/= inEdge) initialEdgeVect
        newEdgeVect = cleanupEdges V.++ newEdges
        newTree = (newVertexVect, newEdgeVect)
        -- add new costs from added vertex to each reamaining leaf
        augmentedDistMatrix = getNewDistMatrix distMatrix addCost eVertLeafDist uVertLeafDist eVertex uVertex leaf
        newDelta = addCost + eVertLeafDist + uVertLeafDist - inWeight
    in
    if (newDelta < inDelta) then (newDelta, newTree, augmentedDistMatrix)
    else addToEdgeSwapRecurse inDelta distMatrix leaf initialTree newLeafIndex (V.tail inEdgeVect)

-- | getVectorAllVectorPairs takes two vectors and creates a vector of avector of two elements each for each
-- pairwise combinatrion of elements
getVectorAllVectorPairs :: V.Vector a -> V.Vector a -> M.Matrix a
getVectorAllVectorPairs firstVect secondVect =
  if V.null firstVect then V.empty
  else
    let firstElement = V.head firstVect
        firstPairs = V.map (V.cons firstElement) $ V.map V.singleton secondVect
    in
    firstPairs V.++ getVectorAllVectorPairs (V.tail firstVect) secondVect

-- | createVectorEdgePairs creates teh Vector of Vectors of edges (2 in each case) to connect 
-- if SPR then takes the initial (previous e edge) and pairs with all in u edge vect
-- if TBR then all conbinations of pairs
createVectorEdgePairs :: String -> V.Vector Edge -> V.Vector Edge -> V.Vector Edge -> V.Vector (V.Vector Edge)
createVectorEdgePairs pairSet previousEdges eEdgeVect uEdgeVect
  | pairSet == "spr" =
    let eEdgePrev = V.head previousEdges
    in
    V.map (V.cons eEdgePrev) $ V.map V.singleton uEdgeVect
  | pairSet == "tbr" = getVectorAllVectorPairs eEdgeVect uEdgeVect
  | otherwise = error ("Pair set option " ++ pairSet ++ " not implemented")


-- | addEdgeToSplitRecurse like addToEdgeSplit but recursiblye yeilds a single best tree
addEdgeToSplitRecurse :: V.Vector Edge -> V.Vector Edge -> Edge -> Edge -> M.Matrix Double -> V.Vector (V.Vector Edge) -> (Double, Tree, M.Matrix Double) -> (Double, Tree, M.Matrix Double)
addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix edgesToConnectVect origTriple@(inDelta, _, _) =
  if V.null edgesToConnectVect then origTriple
  else 
    let edgesToConnect = V.head edgesToConnectVect
    in
    if V.null eEdges &&  V.null uEdges then error "Empty e/u Edge vectors in addEdgeToSplit"
    else if V.null edgesToConnect then error "Empty eEdge vector in addEdgeToSplit"
    else if fst3 (V.head eEdges) == (-1) then -- adding eOTU, V.last since the (-1,-1,0) edge should always be first
      let (newDelta, newTree, newMatrix) = addToEdgeSwap distMatrix (fst3 eTerminal) (V.empty, uEdges) (M.rows distMatrix) (V.last edgesToConnect)
      in 
      if (newDelta < inDelta) then (newDelta, newTree, newMatrix)
      else addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix (V.tail edgesToConnectVect) origTriple
    else if fst3 (V.head uEdges) == (-1) then -- adding uOTU
      let (newDelta, newTree, newMatrix) = addToEdgeSwap distMatrix (fst3 uTerminal) (V.empty, eEdges) (M.rows distMatrix) (V.last edgesToConnect)
      in 
      if (newDelta < inDelta) then (newDelta, newTree, newMatrix)
      else addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix (V.tail edgesToConnectVect) origTriple
    else -- both internal edges
      let (newDelta, newTree, newMatrix) = connectEdges distMatrix eEdges uEdges edgesToConnect
      in 
      if (newDelta < inDelta) then (newDelta, newTree, newMatrix)
      else addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix (V.tail edgesToConnectVect) origTriple


-- | doSPRTBR takes split tree and rejoins by creating edges from a single one of the first edge set to each of the mebers of the second
-- important that e and u edges come in correct order and previous edges are e and u sets in order as well
doSPRTBR :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> [TreeWithData]
doSPRTBR rejoinType curBestCost leafNames outGroup  split =
  -- trace ("In doSPRTBR") (
  let (eEdgeVect, uEdgeVect, delta, previousEdges, distMatrix) = split
  in
  if V.null eEdgeVect || V.null uEdgeVect then error "Empty edge vectors in doSPRTBR"
  else
      if fst3 (V.head eEdgeVect) == snd3 (V.head eEdgeVect) then -- if an OTU call readdTerminal
        let newLeafIndex = M.rows distMatrix
            -- keep whole uEdgeVect so recheck input tree edges
            additionList = V.map (addToEdgeSwap distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex) uEdgeVect
            minAdditionCost = V.minimum (V.map fst3 additionList)
        in
        if minAdditionCost > delta then []
        else sieveTrees delta curBestCost additionList leafNames outGroup []
      else -- internal edge or edge with two OTUs as vertices
          -- check to make sure edges are where they should be can remove later
          -- fix e edge and join to each u edge for SPR (n^2)
          let edgesToConnect = createVectorEdgePairs rejoinType previousEdges eEdgeVect uEdgeVect
              -- e and u terminal should not be used here since OTUs are shortcircuited above
              eTerminal = (-1,-1,0)
              uTerminal = (-1,-1,0)
              additionList = V.map (addEdgeToSplit eEdgeVect uEdgeVect eTerminal uTerminal distMatrix) edgesToConnect
              minAdditionCost = V.minimum (V.map fst3 additionList)
          in
          if minAdditionCost > delta then []
          else sieveTrees delta curBestCost additionList leafNames outGroup []


-- | doSPRTBRSteep like doSPRTBR but only saves a sinlge (and better) tree
doSPRTBRSteep :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData
doSPRTBRSteep rejoinType curBestCost leafNames outGroup split origTree@(_, inTree, _, inMatrix) =
  -- trace ("In doSPRTBR") (
  let (eEdgeVect, uEdgeVect, delta, previousEdges, distMatrix) = split
  in
  if V.null eEdgeVect || V.null uEdgeVect then error "Empty edge vectors in doSPRTBR"
  else
      if fst3 (V.head eEdgeVect) == snd3 (V.head eEdgeVect) then -- if an OTU call readdTerminal
        let newLeafIndex = M.rows distMatrix
            -- keep whole uEdgeVect so recheck input tree edges
            (newDelta, newTree, newMatrix) = addToEdgeSwapRecurse delta distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex uEdgeVect
            newCost = curBestCost - delta + newDelta
            newickTree = convertToNewick leafNames outGroup newTree ++ "[" ++ show newCost ++ "]" ++ ";"
        in
        if (newCost < curBestCost) then (newickTree, newTree, newCost, newMatrix)
        else origTree
      else -- internal edge or edge with two OTUs as vertices
          -- check to make sure edges are where they should be can remove later
          -- fix e edge and join to each u edge for SPR (n^2)
          let edgesToConnect = createVectorEdgePairs rejoinType previousEdges eEdgeVect uEdgeVect
              -- e and u terminal should not be used here since OTUs are shortcircuited above
              eTerminal = (-1,-1,0)
              uTerminal = (-1,-1,0)
              (newDelta, newTree, newMatrix) = addEdgeToSplitRecurse eEdgeVect uEdgeVect eTerminal uTerminal distMatrix edgesToConnect (delta, inTree, inMatrix)
              newCost = curBestCost - delta + newDelta
              newickTree = convertToNewick leafNames outGroup newTree ++ "[" ++ show newCost ++ "]" ++ ";"
          in
          if (newCost < curBestCost) then trace ("->" ++ show newCost) (newickTree, newTree, newCost, newMatrix)
          else origTree


-- | reAddTerminalsSteep like readdTerminals but only returns one tree keeping better 
reAddTerminalsSteep :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData
reAddTerminalsSteep rejoinType curBestCost leafNames outGroup split origTree =
  if rejoinType /= "otu" then error ("Incorrect swap function in reAddTerminals: " ++ rejoinType)
  else
    let (eEdgeVect, uEdgeVect, delta, _, distMatrix) = split
        nOTUs = V.length leafNames
    in
    if (V.length eEdgeVect > 1) || ((fst3 (V.head eEdgeVect) /= snd3 (V.head eEdgeVect)) && (fst3 (V.head eEdgeVect) < nOTUs) && (snd3 (V.head eEdgeVect) < nOTUs)) then origTree
    else if M.rows distMatrix /= ((2 * nOTUs) - 3) then error ("Dist Matrix incorrect size " ++ show (M.dim distMatrix) ++ " should be " ++ show ((2 * nOTUs) - 3, (2 * nOTUs) - 3))
    else
      let newLeafIndex = M.rows distMatrix
      -- take tail of uEdgeVect so not regerate input tree
          (newDelta, newTree, newMatrix) = addToEdgeSwapRecurse delta distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex uEdgeVect -- (V.tail uEdgeVect) -- tail so not hit original tree, leave all to reestimate if necesary
          newCost = curBestCost - delta + newDelta
          newickTree = convertToNewick leafNames outGroup newTree ++ "[" ++ show newCost ++ "]" ++ ";"
      in
      if (newCost < curBestCost) then trace ("->" ++ show newCost) (newickTree, newTree, newCost, newMatrix)
      else origTree



-- | filterNewTreesOnCost returns list of all unique new best cost trees from list
-- assumes curBestCost = cost of sabed trees 
filterNewTreesOnCost :: Double -> [TreeWithData] -> [TreeWithData] -> [TreeWithData]
filterNewTreesOnCost curBestCost firstTreeList savedTrees =
  if null firstTreeList then savedTrees
  else
      let firstTree = head firstTreeList
          (_, _, firstCost, _) = firstTree
      in
      if firstCost < curBestCost then filterNewTreesOnCost firstCost (tail firstTreeList) [firstTree]
      else if firstCost > curBestCost then filterNewTreesOnCost curBestCost (tail firstTreeList) savedTrees
      else
          let uniqueTree = filterNewTrees savedTrees firstTree
          in
          if isNothing uniqueTree then filterNewTreesOnCost curBestCost (tail firstTreeList) savedTrees
          else filterNewTreesOnCost curBestCost (tail firstTreeList) (fromJust uniqueTree : savedTrees )

-- | filterNewTrees takes the first tree and checks if in the second list 
filterNewTrees :: [TreeWithData] -> TreeWithData -> Maybe TreeWithData
filterNewTrees secondTreeList firstTree =
  if null secondTreeList then Just firstTree
  else
    let (firstNewick, _, _, _) = firstTree
        (secondNewick, _, _, _) = head secondTreeList
    in
    if firstNewick == secondNewick then Nothing
    else filterNewTrees (tail secondTreeList) firstTree

-- | getSaveNumber returns number ot save or Infty
getSaveNumber :: String -> Int
getSaveNumber inString =
  if length inString == 4 then maxBound :: Int
  else (read $ drop 5 inString) :: Int

-- | splitJoin does both split and rejoin operations in a fashion that if a better (shorter) tree is found is shortcircuits and 
-- begins again on the new tree, else proceeds untill all splits and joins are completed, but only on a single tree
splitJoin :: (String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData) -> String -> V.Vector String -> Int -> V.Vector Edge -> TreeWithData -> TreeWithData
splitJoin swapFunction refineType leafNames outGroup edgeVect curTreeWithData@(_, curTree, curTreeCost, curTreeMatrix) = 
  if V.null edgeVect then curTreeWithData -- All splits tested, nothing better found
  else 
    let firstEdge = V.head edgeVect
        firstSplit = splitTree curTreeMatrix curTree curTreeCost firstEdge
        !firstTree@(_, firstNewTree, firstTreeCost, _) = swapFunction refineType curTreeCost leafNames outGroup firstSplit curTreeWithData
    in
    if firstTreeCost < curTreeCost then splitJoin swapFunction refineType leafNames outGroup (snd firstNewTree) firstTree
    else splitJoin swapFunction refineType leafNames outGroup (V.tail edgeVect) curTreeWithData


-- | splitJoinWrapper wraps around splitJoin to allow parallel execution
-- the reason is to allow the consumption of "edgeVect" recursively within the same tree
splitJoinWrapper :: (String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData) -> String -> V.Vector String -> Int -> TreeWithData -> TreeWithData
splitJoinWrapper swapFunction refineType leafNames outGroup curTreeWithData@(_, curTree, _, _) =
    let edgeVect = snd curTree
    in
    splitJoin swapFunction refineType leafNames outGroup edgeVect curTreeWithData 

-- | getGeneralSwapSteepestOne performs refinement as in getGeneralSwap but saves on a single tree (per split/swap) and 
-- immediately accepts a Better (shorter) tree and resumes the search on that new tree
-- relies heavily on laziness of splitTree so not parallel at this level
getGeneralSwapSteepestOne :: String -> (String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData) -> V.Vector String -> Int -> [TreeWithData] -> [TreeWithData] -> [TreeWithData]
getGeneralSwapSteepestOne refineType swapFunction leafNames outGroup inTreeList savedTrees =
  if null inTreeList then savedTrees
  else
      trace ("In "++ refineType ++ " Swap (steepest) with " ++ show (length inTreeList) ++ " trees with minimum length " ++ show (minimum $ fmap thd4 inTreeList)) (
      let steepTreeList = seqParMap myStrategy (splitJoinWrapper swapFunction refineType leafNames outGroup) inTreeList
          steepCost = minimum $ fmap thd4 steepTreeList 
      in
      --this to maintina the trajectories untill final swap--otherwise could converge down to single tree prematurely
      keepTrees steepTreeList "unique" "first" steepCost
      {-
      -- saving equal here so can be sent on to full equal tree refine later if nothing better is found
      if steepCost < overallBestCost then getGeneralSwapSteepestOne refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) [steepTree]
      else if steepCost == overallBestCost then getGeneralSwapSteepestOne refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) (steepTree : savedTrees)
      else getGeneralSwapSteepestOne refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees
      -}
      )

-- | getGeneralSwap performs a "re-add" of terminal identical to wagner build addition to available edges
-- performed on all splits recursively until no more better/equal cost trees found
-- this won't work to save "all" just unique and best and unique of best
-- add "steep-est" descent
getGeneralSwap :: String -> (String -> Double -> V.Vector String -> Int -> SplitTreeData -> [TreeWithData]) -> String -> String -> V.Vector String -> Int -> [TreeWithData] -> [TreeWithData] -> [TreeWithData]
getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup inTreeList savedTrees =
  let maxNumSave = getSaveNumber saveMethod
  in
  if null inTreeList then savedTrees
  else
      trace ("In "++ refineType ++ " Swap with " ++ show (length inTreeList) ++ " trees with minimum length " ++ show (minimum $ fmap thd4 inTreeList)) (
      let curFullTree = head inTreeList
          overallBestCost = minimum $ fmap thd4 savedTrees
          (_, curTree, curTreeCost, curTreeMatrix) = curFullTree
          -- parallelize here 
          -- splitTreeList = fmap (splitTree curTreeMatrix curTree curTreeCost) (snd curTree) -- `using` parListChunk chunkSize rdeepseq
          splitTreeList = seqParMap myStrategy (splitTree curTreeMatrix curTree curTreeCost) (snd curTree)
          -- (chunkSize, _) = quotRem (length splitTreeList) getNumThreads
          -- firstTreeList = fmap (swapFunction refineType curTreeCost leafNames outGroup nOTUs) splitTreeList  `using` parListChunk chunkSize rdeepseq
          firstTreeList = seqParMap myStrategy (swapFunction refineType curTreeCost leafNames outGroup) splitTreeList
          -- firstTreeList = V.map (reAddTerminals curBestCost leafNames outGroup nOTUs) splitTreeList
          firstTreeList' = filterNewTreesOnCost overallBestCost  (curFullTree : concat (V.toList firstTreeList)) savedTrees -- keepTrees (concat $ V.toList firstTreeList) saveMethod overallBestCost
      in
      -- trace (show splitTreeList) (
      -- Work around for negative NT.infinity tree costs (could be dst matrix issue)
      if NT.isInfinite curTreeCost || null firstTreeList' then getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees else (
   let (_, _, costOfFoundTrees, _) = head firstTreeList'
   in
   --workaround for negatrive NT.infinity trees
   if NT.isInfinite costOfFoundTrees then getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees
   else if costOfFoundTrees < overallBestCost then
     let uniqueTreesToAdd = fmap fromJust $ filter (/= Nothing ) $ fmap (filterNewTrees inTreeList) firstTreeList'
         treesToSwap = keepTrees (tail inTreeList ++ uniqueTreesToAdd) saveMethod keepMethod costOfFoundTrees
     in
     getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup treesToSwap (take maxNumSave firstTreeList')
   else if costOfFoundTrees == overallBestCost then
     if length savedTrees >= maxNumSave then getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees
     else
      let uniqueTreesToAdd = fmap fromJust $ filter (/= Nothing ) $ fmap (filterNewTrees inTreeList) firstTreeList'
          treesToSwap = keepTrees (tail inTreeList ++ uniqueTreesToAdd) saveMethod keepMethod costOfFoundTrees
      in
      getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup treesToSwap (take maxNumSave $ savedTrees ++ firstTreeList')
   else getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees)
        ) -- )

-- | performRefinement takes input trees in TRE and Newick format and performs differnt forma of tree refinement
-- at present just OTU (remove leaves and re-add), SPR and TBR
performRefinement :: String -> String -> String -> V.Vector String -> Int -> TreeWithData -> [TreeWithData]
performRefinement refinement saveMethod keepMethod leafNames outGroup inTree
  | refinement == "none" = [inTree]
  | refinement == "otu" =
    let !newTrees = getGeneralSwapSteepestOne "otu" reAddTerminalsSteep leafNames outGroup [inTree] [([],(V.empty,V.empty), NT.infinity, M.empty)]
        !newTrees' = getGeneralSwap "otu" reAddTerminals saveMethod keepMethod leafNames outGroup newTrees [([],(V.empty,V.empty), NT.infinity, M.empty)]
    in
    if not (null newTrees') then newTrees'
    else
      trace "OTU swap did not find any new trees"
      [inTree]
  | refinement == "spr" =
    let !newTrees = getGeneralSwapSteepestOne "otu" reAddTerminalsSteep leafNames outGroup [inTree] [([],(V.empty,V.empty), NT.infinity, M.empty)]
        !newTrees' = getGeneralSwapSteepestOne "spr" doSPRTBRSteep leafNames outGroup newTrees [([],(V.empty,V.empty), NT.infinity, M.empty)]
        !newTrees'' = getGeneralSwap "spr" doSPRTBR saveMethod keepMethod leafNames outGroup newTrees' [([],(V.empty,V.empty), NT.infinity, M.empty)]
    in
    if not (null newTrees'') then newTrees''
    else
      trace "SPR swap did not find any new trees"
      [inTree]
  | refinement == "tbr" =
    let !newTrees = getGeneralSwapSteepestOne "otu" reAddTerminalsSteep leafNames outGroup [inTree] [([],(V.empty,V.empty), NT.infinity, M.empty)]
        !newTrees' = getGeneralSwapSteepestOne "spr" doSPRTBRSteep leafNames outGroup newTrees [([],(V.empty,V.empty), NT.infinity, M.empty)]
        !newTrees'' = getGeneralSwapSteepestOne "tbr" doSPRTBRSteep leafNames outGroup newTrees' [([],(V.empty,V.empty), NT.infinity, M.empty)]
        !newTrees''' = getGeneralSwap "tbr" doSPRTBR saveMethod keepMethod leafNames outGroup newTrees'' [([],(V.empty,V.empty), NT.infinity, M.empty)]
    in
    if not (null newTrees''') then newTrees'''
    else
      trace "TBR swap did not find any new trees"
      [inTree]
  | otherwise = error ("Unrecognized refinement method: " ++ refinement)

-- | makeVertexNames takes vertgex indices and returns leaf name if < nOTUs and "HTU" ++ show Index
-- if not 
makeVertexNames :: [Vertex] -> Int -> V.Vector String -> [String]
makeVertexNames vertList nOTUs leafNames =
  if null vertList then []
  else
      let firstVert = head vertList
      in
      if firstVert < nOTUs then (leafNames V.! firstVert) : makeVertexNames (tail vertList) nOTUs leafNames
      else ("HTU" ++ show firstVert) : makeVertexNames (tail vertList) nOTUs leafNames

-- | directSingleEdge takes an Int and makes that 'e' and otehr vertex as 'u' in edge (e->u)
directSingleEdge :: Int -> Edge -> Edge
directSingleEdge index (a,b,w)
  | a == index = (a,b,w)
  | b == index = (b,a,w)
  | otherwise = error ("Index " ++ show index ++ " doesn't match edge " ++ show (a,b,w) ++ " in directSingleEdge")

-- | getChildEdges returns the two edges that are childre of a vertex
getChildEdges :: Int -> Int -> V.Vector Edge -> V.Vector Edge
getChildEdges vertIndex nLeaves inEdgeVect
  | V.null inEdgeVect = V.empty
  | vertIndex < nLeaves = error ("Looking for child of leaf " ++ show (vertIndex, nLeaves))
  | otherwise =
    let (a,b,w) = V.head inEdgeVect
    in
    if (a == vertIndex) || (b == vertIndex) then V.cons (a,b,w) (getChildEdges vertIndex nLeaves (V.tail inEdgeVect)) else getChildEdges vertIndex nLeaves (V.tail inEdgeVect)

-- | directexEdges takes a vector of edges and outgrop index and directs the edges (parent -> child vertices) based on that
directEdges :: Int -> Int -> Bool -> V.Vector Edge -> V.Vector Edge
directEdges vertIndex nLeaves isFirst inEdgeVect
  | V.null inEdgeVect = V.empty
  | isFirst = --to find out group edge order larger to smaller will have outgroup index second
    let outgroupEdge = getEdgeRoot vertIndex inEdgeVect
        remainingEdgeVect = subtractVector (V.singleton outgroupEdge) inEdgeVect
        (a,b,w) = orderEdge outgroupEdge
    in
    V.cons (a,b,w) (directEdges a nLeaves False remainingEdgeVect)
  | vertIndex < nLeaves = V.empty
  | otherwise = -- not outgroup but regular node, get two child edges
    let descdendantEdges = getChildEdges vertIndex nLeaves inEdgeVect
        remainingEdgeVect = subtractVector descdendantEdges inEdgeVect
        newDescEdges = V.map (directSingleEdge vertIndex) descdendantEdges
    in
    if V.length newDescEdges /= 2 then error ("There should be 2 child edges for index " ++ show vertIndex ++ " and there are(is) " ++ show (V.length newDescEdges) ++ " " ++ show newDescEdges)
    else
        let (_, bf, _) = V.head newDescEdges
            (_, bs, _) = V.last newDescEdges
            firstSubEdges = directEdges bf nLeaves False remainingEdgeVect
            remainingEdgeVect' = subtractVector firstSubEdges remainingEdgeVect
            secondSubEdges = directEdges bs nLeaves False remainingEdgeVect'
        in
        (newDescEdges V.++ (firstSubEdges V.++ secondSubEdges))


-- | convertToGraph takes Vertex of names and a tree and return inductive Graph format
convertToDirectedGraph :: V.Vector String -> Int -> Tree -> P.Gr String Double
convertToDirectedGraph leafList outgroupIndex inTree =
  let (_, edgeVect) = inTree
      nOTUs = length leafList
      -- should be stright 0->n-1 but in case some vertex number if missing
      vertexList = sort $ Set.toList $ getVertexSet edgeVect
      vertexNames = makeVertexNames vertexList nOTUs leafList
      labelledVertexList = Data.List.zip vertexList vertexNames
      edgeList = V.toList $ directEdges outgroupIndex nOTUs True edgeVect
  in
  G.mkGraph labelledVertexList edgeList

-- | showGraph a semi-formatted show for Graphs
-- showGraph :: P.Gr String Double ->  String
showGraph :: (Show a, Show b) => P.Gr a b -> String
showGraph inGraph =
  if G.isEmpty inGraph then "Empty Graph"
  else
      let nodeString = show $ G.labNodes inGraph
          edgeString  = show $ G.labEdges inGraph
      in
      ("Nodes:" ++ nodeString ++ "\n" ++ "Edges: " ++ edgeString)

-- | writeFiles takes a stub and list of strings
-- and writes to files usinf stub ++ number as naming convention
writeFiles :: String -> String -> Int -> [String] -> IO ()
writeFiles stub typeString number fileStuffList =
  if null fileStuffList then hPutStrLn stderr ("Wrote " ++ show number ++ " " ++ stub ++ ".X." ++ typeString ++ " files")
  else do
    writeFile (stub ++ "." ++ show number ++ "." ++ typeString) (head fileStuffList)
    writeFiles stub typeString (number + 1) (tail fileStuffList)

-- | removeCommenet removes all after and including double dash "--"
removeComments :: String -> String
removeComments inLine
  | null inLine = []
  | length inLine == 1 = inLine
  | otherwise =
    let firstElem = head inLine
        secondElem  = inLine !! 1
    in
    if not ((firstElem == '-') && (secondElem  == '-'))  then firstElem : removeComments (tail inLine)
    else []

-- | cleanUpParamFile removes comment lines and space etc from input parameter file
-- retuns clean lines as lists of string
cleanUpParamFile :: String -> [String]
cleanUpParamFile inFile =
  if null inFile then error "Empty input file to clean"
  else
      let inLines = lines inFile
          nonCommentLines = filter (/= "") $ fmap (filter (/= ' ') . removeComments) inLines
      in
      nonCommentLines

-- | getOption scanns input list of command strings for the string input and returns that line
-- for later parsing, converts all to lower caser so case invariant later
-- uses dataFiles String as default for stub if not specified
getOption :: String -> String -> [String] -> String
getOption optionString dataFileString commandLines =
  if null commandLines then
    -- return defaults or error
   case optionString of
        "input"             -> error "No input data file specified"
        "stub"              -> dataFileString
        "output"            -> dataFileString ++ ".tre"
        "firstPairChoice"   -> "closest"
        "outgroup"          -> []
        "additionSequence"  -> "best"
        "refinement"        -> "none"
        "buildSet"          -> "best"
        "outputSet"         -> "best"
        "keepSet"           -> "first"
        "excludedTaxa"      -> []
        _                   -> error ("Option " ++ optionString ++ " not specified and has no default")
  else -- case invariant
      let firstLine = head commandLines
          parameterString = tail $ dropWhile (/= ':') firstLine
          parameterStringLC = T.unpack $ T.toLower $ T.pack parameterString
          commandString = T.toLower $ T.pack $ takeWhile (/= ':') firstLine
          inOption = T.toLower $ T.pack optionString
      in
      if inOption == commandString then
          if optionString `elem` ["input","output","stub","outgroup","excludedTaxa"] then parameterString
          else  parameterStringLC
      else getOption optionString dataFileString (tail commandLines)

-- | processParamFile takes input Wagner script file and returns run parameters
-- as a list if string in order
-- input Data File, firstPairMethod, outgroup, additionSeqeunce, refinement, buildSelect, saveMethod, stub,
-- outputTreeFile, keep method
processParamFile :: String -> [String]
processParamFile fileString =
  if null fileString then error "Empty parameter file"
  else
      let commandLines      = cleanUpParamFile fileString
          inputFile         = filter (/= '"') $ getOption "input" "" commandLines
          firstPair         = getOption "firstPairChoice" inputFile commandLines
          outgroup          = getOption "outgroup" inputFile commandLines
          additionSequence  = getOption "additionSequence" inputFile commandLines
          refinement        = getOption "refinement" inputFile commandLines
          buildSelect       = getOption "buildSet" inputFile commandLines
          saveMethod        = getOption "outputSet" inputFile commandLines
          stub              = filter (/= '"') $ getOption "stub" inputFile commandLines
          outputTreeFile    = getOption "output" inputFile commandLines
          keepMethod        = getOption "keepSet" inputFile commandLines
          excludedTaxa      = getOption "excludedTaxa" inputFile commandLines
      in
      -- check for unset options-> throw error
      if null inputFile then error "No input file specified"

      else [inputFile, firstPair, outgroup, additionSequence, refinement, buildSelect, saveMethod, stub, outputTreeFile, keepMethod,excludedTaxa]

-- |  getDeletedTaxa returns list of taxa to delete from file exludedTaxaFileName
getDeletedTaxa :: String -> IO [String]
getDeletedTaxa exludedTaxaFileName =
  if null exludedTaxaFileName then return []
  else
    do
      fileContents <- readFile $ filter (/= '"') exludedTaxaFileName
      return $ words fileContents

-- makeNewRow takes a list of indices and deltes those elements form a list of Strings
makeNewRow :: [Int] -> [String] -> Int -> [String]
makeNewRow taxonIndexList stringRow columnNum
  | null stringRow = []
  | columnNum `elem` taxonIndexList = makeNewRow taxonIndexList (tail stringRow) (columnNum + 1)
  | otherwise =
    head stringRow : makeNewRow taxonIndexList (tail stringRow) (columnNum + 1)

-- makeNewData take rawData and list of rows and columns to delete
makeNewData :: [Int] -> [[String]] -> Int -> [[String]]
makeNewData taxonIndexList inData rowNum
  | null taxonIndexList = inData
  | null inData = []
  | rowNum `elem` taxonIndexList = makeNewData taxonIndexList (tail inData) (rowNum + 1)
  | otherwise =
    let firstRow = makeNewRow taxonIndexList (head inData) 0
    in
    firstRow : makeNewData taxonIndexList (tail inData) (rowNum + 1)


-- | deleteTaxa deletes taxa from rawData (from csv) and returns new matrixc
-- with top line of taxon names
deleteTaxa :: [String] -> [[String]] -> [[String]]
deleteTaxa toDeleteList inData
  | null inData = error "Empty input data in deleteTaxa (empty data file)"
  | null toDeleteList = inData
  | otherwise =
    let taxonList = head inData
        fullTaxonIndexList = fmap (`elemIndex` taxonList) toDeleteList
        taxonIndexList = fromJust <$> filter (/= Nothing) fullTaxonIndexList
        firstRow = taxonList \\ toDeleteList
    in
    if null taxonIndexList then inData
    else -- delete rows and columns, first row done above (name row)
      firstRow : makeNewData taxonIndexList (tail inData) 0

-- | getOutgroup sets outgroup to String input or 0 (first taxon) by default.
getOutgroup :: String -> V.Vector String -> Int
getOutgroup outString leafNames =
  if null outString then
      trace ("\nDefault Rooting on first leaf : " ++  (leafNames V.! 0))
      0
  else
    let outElem = V.elemIndex outString leafNames
    in
    if isNothing outElem then error ("Outgroup name " ++ outString ++ " not found")
    else
      trace ("\nRooting on leaf: " ++ outString)
      fromJust outElem

-- | main driver
main :: IO ()
main =
  do
    -- Process arguments
    args <- getArgs
    if length args /= 1 then error "Need to specify a single parameter file"
    else hPutStrLn stderr ("Openning parameter file " ++ head args)

    -- Get params from input file
    paramFile <- readFile (head args)
    let paramList = processParamFile paramFile
    let dataFile = head paramList
    let firstPairMethod = paramList !! 1 -- should be closest, furthest, random
    let outgroup = filter (/= '"') $ paramList !! 2
    let addSequence =  paramList !! 3
    let refinement = paramList !! 4
    let buildSelect = paramList !! 5
    let saveMethod = paramList !! 6
    let stub = paramList !! 7
    let outputTreeFile = filter (/= '"')  $ paramList !! 8 -- remove quotes
    let keepMethod =  paramList !! 9
    let exludedTaxaFileName =  paramList !! 10 -- remove quotes

    -- Prelude.mapM_ (hPutStrLn stderr) args
    csvResult <- parseFromFile csvFile dataFile
    let rawDataInit = case csvResult of
                      Left err -> error $ "Error parsing " ++ dataFile ++ " " ++ show err
                      Right result -> result
    let rawData = filter (/=[""]) rawDataInit -- remove crap from csv

    --Delete taxa to be excluded
    leavesToDelete <- getDeletedTaxa exludedTaxaFileName
    hPutStrLn stderr "Deleted taxa "
    Prelude.mapM_ (hPutStrLn stderr) leavesToDelete
    let rawData' = deleteTaxa leavesToDelete rawData
    -- hPutStrLn stderr $ show rawData'

    let leafNames = V.fromList $ head rawData'
    hPutStrLn stderr "\nRemaining  taxa "
    Prelude.mapM_ (hPutStrLn stderr) (fmap show leafNames)

    -- Set outgroup with default to '0'
    let !outElem = getOutgroup outgroup leafNames

    if head addSequence /= 'r' then hPutStrLn stderr ("\nAddition sequence: " ++ addSequence )
    else hPutStrLn stderr ("\nAddition sequence: random with " ++ drop 7 addSequence ++ " replicates" )
    hPutStrLn stderr ("Output set: " ++ saveMethod )

    if V.length leafNames /= length (tail rawData') then error ("Input matrix is not square: " ++ show (V.length leafNames) ++ " leaves and " ++
      show (length $ tail rawData') ++ " rows")
    else hPutStrLn stderr ("Input distance data for " ++ show (V.length leafNames) ++ " leaves")
    let rowLengthCheck = foldl' (&&) True $ fmap ((== V.length leafNames) . length) (tail rawData')
    if not rowLengthCheck then error "Row lengths do not equal leaf number"
    else hPutStrLn stderr "Input matrix is square"

    -- Convert to matrix of Doubles
    let distMatrix = M.fromLists $ fmap (fmap (read :: String -> Double)) (tail rawData')

    -- Callan random shuffle
    let randomAddsToDo = getRandomReps addSequence
    let testLeavesVect = V.fromList [0..(V.length leafNames - 1)]
    -- let (chunkSize, _) = quotRem randomAddsToDo getNumThreads
    shuffledList <- CMP.replicateM randomAddsToDo (shuffleM testLeavesVect) -- `using` parListChunk chunkSize rdeepseq

    let !treeList = doWagnerS leafNames distMatrix firstPairMethod outElem addSequence shuffledList

    -- Filter trees from build
    let filteredTrees = keepTrees treeList buildSelect keepMethod NT.infinity -- modify to keep Tree best as well

    hPutStrLn stderr ("After build, there are " ++ show (length filteredTrees) ++ " saved trees at cost " ++ show (minimum $ fmap thd4 filteredTrees))


    let !refinedTrees = concat $ seqParMap myStrategy (performRefinement refinement saveMethod keepMethod leafNames outElem) filteredTrees

    --finale keep 
    let finalTrees = keepTrees refinedTrees saveMethod keepMethod NT.infinity

    hPutStrLn stderr ("After refinement, there are " ++ show (length finalTrees) ++ " saved trees at cost " ++ show (thd4 $ head finalTrees))

    -- output newick file
    hPutStrLn stderr ("Writing output tree file " ++ outputTreeFile)
    writeFile outputTreeFile (unlines (fmap fst4 finalTrees))

    -- Output "dot" for graphviz
      -- Convert to Graph (fgl) format 
    let graphList = fmap (convertToDirectedGraph leafNames outElem . snd4) finalTrees

      -- Convert to Dot format
    let dotList = fmap (GV.graphToDot GV.quickParams) graphList

      -- Output Dot 
    let dotStringList = fmap ((T.unpack . renderDot) . toDot) dotList


    writeFiles stub "dot" 0 dotStringList
    -- writeFiles (reverse $ tail $ dropWhile (/= '.') (reverse dataFile)) "dot" 0 dotStringList

    hPutStrLn stderr "Done"


