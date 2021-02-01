{- |
Module      :  ParseCommands.hs
Description :  Progam to parseom commands from commandline or file
               input graphviz dot files and newick
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

Todo:
  from input file?
  add edit distance to commands for error

-}

module ParseCommands (processParamString) where

import qualified Data.Text.Lazy         as T
-- import           Debug.Trace
import           Data.Array

-- | editDistance is a naive edit distance between two lists
-- takes two  lists and returns edit distance
--- from  https://wiki.haskell.org/Edit_distance
editDistance :: Eq a => [a] -> [a] -> Int
editDistance xs ys = table ! (m,n)
    where
    (m,n) = (length xs, length ys)
    x     = array (1,m) (zip [1..] xs)
    y     = array (1,n) (zip [1..] ys)

    table :: Array (Int,Int) Int
    table = array bnds [(ij, dist ij) | ij <- range bnds]
    bnds  = ((0,0),(m,n))

    dist (0,j) = j
    dist (i,0) = i
    dist (i,j) = minimum [table ! (i-1,j) + 1, table ! (i,j-1) + 1,
        if x ! i == y ! j then table ! (i-1,j-1) else 1 + table ! (i-1,j-1)]

-- | allowedCommandList list of allowable commands
allowedCommandList :: [String]
allowedCommandList = ["input","stub", "output", "firstpairchoice","outgroup",
                    "additionsequence" ,"refinement" ,"buildset" ,"outputset" ,"keepset","excludedtaxa"]

-- allowedOptionsList :: [String]
-- allowedOptionsList = ["closest","best","first","random"]

-- | getBestMatch compares input to allowable commands and checks if in list and if not outputs
-- closest match
getBestMatch :: (Int, String) -> [String] -> String -> (Int, String)
getBestMatch curBest@(minDist, _) allowedStrings inString =
    if null allowedStrings then curBest
    else
        let candidate =  head allowedStrings
            candidateEditCost = editDistance candidate inString
        in
        if candidateEditCost == 0 then (0, candidate)
        else if candidateEditCost < minDist then getBestMatch (candidateEditCost, candidate) (tail allowedStrings) inString
        else getBestMatch curBest (tail allowedStrings) inString

-- | getCommandErrorString takes list of non zero edits to allowed commands and reurns meaningful error string
getCommandErrorString :: [(Int, String, String)] -> String
getCommandErrorString noMatchList =
    if null noMatchList then ""
    else
        let (_, firstCommand, firstMatch) = head noMatchList
            firstError = "\tBy \'" ++ firstCommand ++ "\' did you mean \'" ++ firstMatch ++ "\'?\n"
        in
        firstError ++ getCommandErrorString (tail noMatchList)


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
      if (T.unpack commandString) `notElem` allowedCommandList then
        let (editCost, bestMatch) = getBestMatch (maxBound :: Int ,"no suggestion") allowedCommandList (T.unpack commandString)
        in
        error ("\n\nError in command specification (case insensitive) unrecognized command:" ++ getCommandErrorString [(editCost, T.unpack commandString, bestMatch)])
      else if inOption == commandString then
          if optionString `elem` ["input","output","stub","outgroup","excludedTaxa"] then parameterString
          else  parameterStringLC
      else getOption optionString dataFileString (tail commandLines)

-- | processParamFile takes input Wagner script file and returns run parameters
-- as a list if string in order
-- input Data File, firstPairMethod, outgroup, additionSeqeunce, refinement, buildSelect, saveMethod, stub,
-- outputTreeFile, keep method
processParamString :: String -> Bool -> [String]
processParamString fileString isFile =
  if null fileString then error "Empty parameter file"
  else
      let commandLines      = if isFile then  cleanUpParamFile fileString else words fileString
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



