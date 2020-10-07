{- |
Module      :  Utilities.hs
Description :  Module with useful functionsfor  distance tree construction methods dWag, Neightbor-Joining, UPGMA, and WPGMA
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

-}

module Utilities where

import qualified Data.Vector                       as V 
import qualified SymMatrix                         as M
import           Types
import           Control.Concurrent
import           Control.Parallel.Strategies
import           System.IO.Unsafe
import           Data.List
import qualified System.Random.Shuffle             as RandS
import qualified System.Random                     as Rand
import qualified Data.Set                          as Set



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

fst4 :: (a,b,c,d) -> a
fst4 (e,_,_,_) = e

snd4 :: (a,b,c,d) -> b
snd4 (_,e,_,_) = e

thd4 :: (a,b,c,d) -> c
thd4 (_,_,e,_) = e

-- Parallel stuff
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

-- | getNumThreads gets number of COncurrent  threads
{-# NOINLINE getNumThreads #-}
getNumThreads :: Int
getNumThreads = unsafePerformIO getNumCapabilities

-- | convertToNewick wrapper to remove double commas
convertToNewick :: V.Vector String -> Int -> Tree -> String
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
    in
    if fstNewick `notElem` fmap fst4 uniqueList then getUniqueTrees (tail inList) (firstTree : uniqueList)
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

  -- | ranList generates random list of positive integers
ranList :: Rand.StdGen -> Int -> Int -> [Int]
ranList sg n maxValue = take n $ Rand.randomRs (0,maxValue) sg

-- | driver function to generate list of positive integers
{-# NOINLINE betterRandomList #-}
betterRandomList :: Int -> Int -> [Int]
betterRandomList n maxValue= unsafePerformIO $ do
    sg <- Rand.getStdGen
    return $ ranList sg n maxValue

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

-- | getVertexSet take a vector of edges and creates the set of vertex numbers
getVertexSet :: V.Vector Edge -> Set.Set Vertex
getVertexSet edgeVect =
  if V.null edgeVect then Set.empty
  else
      let (a,b,_) = V.head edgeVect
          thisSet =  Set.fromList [a,b]
      in
      Set.union thisSet (getVertexSet $ V.tail edgeVect)

-- | getMatrixMinPair takes distMatrix initla pinteger pair and value
-- traverses teh matrix and return minimum distance and index pair
-- if tie takes first
  -- call with (-1, -1, NT.infinity) 0 0
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
-- call with (-1 , -1, 0 :: Double) 0 0
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

