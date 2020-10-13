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
{-# LANGUAGE BangPatterns #-}

module Main where

import           Data.CSV
import           Data.List
import           Data.Maybe
import qualified Data.Vector                       as V hiding (replicateM)
import           Debug.Trace
import           System.Environment
import           System.IO
import           Text.ParserCombinators.Parsec
import qualified Control.Monad.Parallel            as CMP
import qualified Data.Graph.Inductive.Graph        as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import qualified Data.GraphViz                     as GV
import           Data.GraphViz.Printing
import qualified Data.Number.Transfinite           as NT
import qualified Data.Set                          as Set
import qualified Data.Text.Lazy                    as T
import           Immutable.Shuffle
import qualified SymMatrix                         as M
import           DistanceMethods
import           Types
import           Utilities
-- import Control.Monad (replicateM)


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

-- | getRandomReps takes teh build arguemnrt and returns the number of random replicates 0 if not ranomd the numbewr otherwise
getRandomReps :: String -> Int
getRandomReps inString
  | head inString /= 'r' = 0
  | length inString < 8 = error ("Need to specify repicate number after '=' in " ++ inString)
  | otherwise = read (drop 7 inString) :: Int




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

    -- Checks on input distances
    if V.length leafNames /= length (tail rawData') then error ("Input matrix is not square: " ++ show (V.length leafNames) ++ " leaves and " ++
      show (length $ tail rawData') ++ " rows")
    else hPutStrLn stderr ("Input distance data for " ++ show (V.length leafNames) ++ " leaves")
    let rowLengthCheck = foldl' (&&) True $ fmap ((== V.length leafNames) . length) (tail rawData')
    if not rowLengthCheck then error "Row lengths do not equal leaf number"
    else hPutStrLn stderr "Input matrix is square"


    -- Convert to matrix of Doubles
    let distMatrix = M.fromLists $ fmap (fmap (read :: String -> Double)) (tail rawData')

    -- Check that distances are non-negative
    let nonNegative =  M.map (>= 0.0) distMatrix
    let nonNegative' = foldl' (&&) True $ V.map (foldl' (&&) True) nonNegative
    _ <- if nonNegative' then
                hPutStrLn stderr "Input distances non-negative"
            else
                error "Input distance has negative values--must all be >= 0.0"

    -- Callan random shuffle
    let randomAddsToDo = getRandomReps addSequence
    let testLeavesVect = V.fromList [0..(V.length leafNames - 1)]
    -- let (chunkSize, _) = quotRem randomAddsToDo getNumThreads
    
    -- Wagner build and refine 
    shuffledList <- CMP.replicateM randomAddsToDo (shuffleM testLeavesVect) -- `using` parListChunk chunkSize rdeepseq
    let !treeList = if (head addSequence) == 'n' then [neighborJoining leafNames distMatrix outElem]
                    else if (head addSequence) == 'w' then [wPGMA leafNames distMatrix outElem]
                    else doWagnerS leafNames distMatrix firstPairMethod outElem addSequence shuffledList

    -- Filter trees from build
    let !filteredTrees = keepTrees treeList buildSelect keepMethod NT.infinity -- modify to keep Tree best as well

    hPutStrLn stderr ("After build, there are " ++ show (length filteredTrees) ++ " saved trees at cost " ++ show (minimum $ fmap thd4 filteredTrees))


    let !refinedTrees = concat $ seqParMap myStrategy (performWagnerRefinement refinement saveMethod keepMethod leafNames outElem) filteredTrees

    --final keep
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

    hPutStrLn stderr "Done"


