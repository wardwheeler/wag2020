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
    bitvectorts?

  Newick smaller subtree left, bigger right for figtree output
    count bits

  Test keep options and swapping on multiple input trees from build

  precision for showing branch lengths could be an input option ie 0 for not showing them at all.

  Clean up code--remove n^2 bad stuff
    TBR time complexity seems pretty awful--if corrrect

  Change matrix to SymMatrixSeq for better time complexity on matrix updates

  Change Vector type to Sequence for better cons/snoc 

-}

module Main where

import qualified Control.Monad.Parallel        as CMP
import           Data.CSV
import qualified Data.GraphViz                 as GV
import           Data.GraphViz.Printing
import           Data.List
import           Data.Maybe
import qualified Data.Number.Transfinite       as NT
import qualified Data.Text.Lazy                as T
import qualified Data.Vector                   as V hiding (replicateM)
import           Debug.Trace
import           DistanceMethods
import           GeneralUtilities
import           Immutable.Shuffle
import           ParallelUtilities
import           ParseCommands
import qualified SymMatrix                     as M
import           System.Environment
import           System.IO
import           Text.ParserCombinators.Parsec
import           Utilities



-- | writeFiles takes a stub and list of strings
-- and writes to files usinf stub ++ number as naming convention
writeFiles :: String -> String -> Int -> [String] -> IO ()
writeFiles stub typeString number fileStuffList =
  if null fileStuffList then hPutStrLn stderr ("Wrote " ++ show number ++ " " ++ stub ++ ".X." ++ typeString ++ " files")
  else do
    writeFile (stub ++ "." ++ show number ++ "." ++ typeString) (head fileStuffList)
    writeFiles stub typeString (number + 1) (tail fileStuffList)

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
  | null inData = errorWithoutStackTrace "Empty input data in deleteTaxa (empty data file)"
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
    if isNothing outElem then errorWithoutStackTrace ("Outgroup name " ++ outString ++ " not found")
    else
      trace ("\nRooting on leaf: " ++ outString)
      fromJust outElem

-- | getRandomReps takes teh build arguemnrt and returns the number of random replicates 0 if not ranomd the numbewr otherwise
getRandomReps :: String -> Int
getRandomReps inString
  | head inString /= 'r' = 0
  | length inString < 8 = errorWithoutStackTrace ("Need to specify repicate number after '=' in " ++ inString)
  | otherwise = read (drop 7 inString) :: Int


-- | main driver
main :: IO ()
main =
  do
    -- Process arguments
    args <- getArgs
    if null args then errorWithoutStackTrace "Need to specify a single parameter file or commandline options"
    else if length args == 1 then hPutStrLn stderr ("Openning parameter file " ++ head args)
    else hPutStrLn stderr ("Processing commandline options: " ++ unlines args)

    -- Get params from input file
    paramFile <- if length args == 1 then readFile (head args) else return ""
    let paramList = if length args == 1 then processParamString paramFile True else processParamString (unwords args) False
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
                      Left err -> errorWithoutStackTrace $ "Error parsing " ++ dataFile ++ " " ++ show err
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
    if V.length leafNames /= length (tail rawData') then errorWithoutStackTrace ("Input matrix is not square: " ++ show (V.length leafNames) ++ " leaves and " ++
      show (length $ tail rawData') ++ " rows")
    else hPutStrLn stderr ("Input distance data for " ++ show (V.length leafNames) ++ " leaves")
    let rowLengthCheck = foldl' (&&) True $ fmap ((== V.length leafNames) . length) (tail rawData')
    if not rowLengthCheck then errorWithoutStackTrace "Row lengths do not equal leaf number"
    else hPutStrLn stderr "Input matrix is square"


    -- Convert to matrix of Doubles
    let distMatrix = M.fromLists $ fmap (fmap (read :: String -> Double)) (tail rawData')

    -- Check that distances are non-negative
    let nonNegative =  M.map (>= 0.0) distMatrix
    let nonNegative' = foldl' (&&) True $ V.map (foldl' (&&) True) nonNegative
    _ <- if nonNegative' then
                hPutStrLn stderr "Input distances non-negative"
            else
                errorWithoutStackTrace "Input distance has negative values--must all be >= 0.0"

    -- Callan random shuffle
    let randomAddsToDo = getRandomReps addSequence
    let testLeavesVect = V.fromList [0..(V.length leafNames - 1)]
    -- let (chunkSize, _) = quotRem randomAddsToDo getNumThreads

    -- Wagner build and refine
    shuffledList <- CMP.replicateM randomAddsToDo (shuffleM testLeavesVect) -- `using` parListChunk chunkSize rdeepseq
    let !treeList
          | head addSequence == 'n' = [neighborJoining leafNames distMatrix outElem]
          | head addSequence == 'w' = [wPGMA leafNames distMatrix outElem]
          | otherwise = doWagnerS leafNames distMatrix firstPairMethod outElem addSequence shuffledList

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


