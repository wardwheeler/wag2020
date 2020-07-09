{- |
Module      :  SymMatrix.hs 
Description :  Progam to manipulate square symmetric lower diagonal matrices with diagnonal values
                as if they were normal matrices
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

module SymMatrix (empty, dim, fromLists, Matrix, 
                   SymMatrix.null, cols, rows,
                   (!), toLists, updateMatrix,
                   addMatrixRows) where

-- import Debug.Trace
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Sort as S

-- | functions for triples
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

-- | Matrix type as Vector of Vectors
type Matrix a = V.Vector (V.Vector a)

-- | empty matrix value from Vector
empty :: Matrix a
empty = G.empty

-- | dim returns dimmension (rows = cols)
dim :: Matrix a -> (Int, Int)
dim inMatrix = (V.length inMatrix, V.length inMatrix)

-- | rows returns number rows in Matrix (rows = cols)
rows :: Matrix a -> Int
rows inM = fst $ dim inM

-- | cols returns number cols in Matrix (rows = cols)
cols :: Matrix a -> Int
cols inM = fst $ dim inM

-- | null returns True of row number is 0
null :: (Eq a) => Matrix a -> Bool
null inMatrix = inMatrix == empty

-- | fromLists takes list of list of a and returns lower diagnoal (with diagonal)
-- matrix as Matrix (Vector of Vectors)
fromLists :: (Eq a, Show a) => [[a]] -> Matrix a
fromLists inListList =
    if L.null inListList then empty
    else 
        let initialSquare = V.map V.fromList $ V.fromList inListList
            colsH = V.length $ V.head initialSquare
            rowsH = V.length initialSquare
        in
        if (colsH /= rowsH) then error ("Input matrix is not square " ++ show inListList)
        else 
            let indexPairs = cartProd [0..(rowsH - 1)] [0..(rowsH - 1)]
                sym = checkSymmetry initialSquare indexPairs
            in
            if not sym then error "Input matrix not symmetrical"
            else makeLowerDiag initialSquare 0 rowsH

-- | toLists takes a Matrix and returns a list of lists (not all same length)
toLists :: (Eq a, Show a) => Matrix a -> [[a]]
toLists inM =
    if SymMatrix.null inM then []
    else 
        V.toList $ V.map V.toList inM

-- | indexing lower diag matrix
(!) :: Matrix a -> (Int, Int) -> a 
(!) inM (iIndex,jIndex) =
    if iIndex > jIndex then 
        (inM V.! iIndex) V.! jIndex
    else
        (inM V.! jIndex) V.! iIndex

-- | makeLowerDiag take a Vector of Vetors (Matrix) and retusn a lower diagonal matrix
-- including diagonal
makeLowerDiag :: (Eq a) => Matrix a -> Int -> Int -> Matrix a
makeLowerDiag inM row numRows=
    if SymMatrix.null inM then error "Input matrix is empty in makeLowerDiag"
    else if (row == numRows) then V.empty
    else 
        let origRow = inM V.! row
            newRow = V.take (row + 1) origRow
        in
        V.cons newRow  (makeLowerDiag inM (row + 1) numRows)

-- | cartesian product of two lists 
cartProd :: [a] -> [a] -> [(a,a)]
cartProd xs ys = [(x,y) | x <- xs, y <- ys]

-- | checkSymmetry takes a Vector of VEctors of a 
-- and checks for symmetry
checkSymmetry :: (Eq a, Show a) => Matrix a -> [(Int, Int)] -> Bool
checkSymmetry inVV pairList =
    if SymMatrix.null inVV then error ("Input matrix is empty")
    else if L.null pairList then True
    else
        let (iIndex, jIndex) = head pairList
            firstCheck = ((inVV V.! iIndex) V.! jIndex) == ((inVV V.! jIndex) V.! iIndex)
        in
        if firstCheck then checkSymmetry inVV (tail pairList)
        else error ("Matrix is not symmetrical:" ++ (show (iIndex, jIndex)) ++ "=>" ++ (show ((inVV V.! iIndex) V.! jIndex)) ++ " /= " ++ (show ((inVV V.! jIndex) V.! iIndex)))

-- | addMatrixRows add rows to existing matrix to extend Matrix dimention
addMatrixRows :: (Eq a) => Matrix a -> Matrix a -> Matrix a
addMatrixRows inM newRows =
    if SymMatrix.null newRows then inM
    else 
        inM V.++ newRows

-- | reIndexTriple taske (i,j,k) and returns (max i j, min i j, k)
reIndexTriple :: (Ord a) => (a, a, b) -> (a, a, b)
reIndexTriple trip@(iIndex, jIndex, value) =
    if iIndex > jIndex then trip
    else (jIndex, iIndex, value)

-- | updateMatrix takes a list of triples and update matrix
-- update all at once checking for bounds
-- could naively do each triple in turn, but would be alot of copying
updateMatrix :: (Eq a, Show a) => Matrix a -> [(Int, Int, a)] -> Matrix a
updateMatrix inM modList =
    if L.null modList then inM
    else 
        let orderedTripleList = S.uniqueSortOn fst3 $ fmap reIndexTriple modList
            minRow = fst3 $ head orderedTripleList
            maxRow = fst3 $ last orderedTripleList
        in
        if maxRow >= rows inM then error ("Update matrix out of bounds, row = " ++ (show $ rows inM) ++ " and trying to update row " ++ (show maxRow))
        else
            let firstPart = V.unsafeTake minRow inM
                restPart  = V.unsafeDrop minRow inM
                modifiedRemainder = updateRows restPart orderedTripleList minRow 
            in
            addMatrixRows firstPart modifiedRemainder

-- | updateRows takes the section of the matrix containing rows that wil be modified
-- (some not) and modifes or copies rows and rerns a Matrix (vector of roow vectors)
-- doesn't do "smart" update of rows, but recursively till no more for that row
updateRows :: (Show a, Eq a) => Matrix a -> [(Int, Int, a)] -> Int -> Matrix a
updateRows inM tripList currentRow =
    if L.null tripList then inM
    else if SymMatrix.null inM then error ("Matrix is empty and there are modifications that remain: " ++ (show $ tripList))
    else 
        let (rowIndex, columnIndex, value) = L.head tripList
            firstOrigRow = V.head inM
        in
        if currentRow /= rowIndex then firstOrigRow `V.cons` (updateRows (V.tail inM) tripList (currentRow + 1))
        else -- account for multiple modifications to same row
            let (newRow, newTripList) = modifyRow firstOrigRow columnIndex value currentRow (L.tail tripList)
            in
            newRow `V.cons` (updateRows (V.tail inM) newTripList (currentRow + 1))

-- | modifyRow takes an initial modification (column and value) and then checks to see if there are more modifications in that
-- row (rowNumber) in the remainder of the list of modifications, returning the new row and mod list as a pair
modifyRow :: V.Vector a -> Int -> a -> Int -> [(Int, Int, a)] -> (V.Vector a, [(Int, Int, a)])
modifyRow inRow colIndex value rowNumber modList =
    if colIndex >= (V.length inRow) then error ("Column to modify is outside length of row " ++ (show $ (rowNumber, colIndex)))
    else 
        let firstPart = V.unsafeTake colIndex inRow
            remainderPart = V.unsafeDrop (colIndex + 1) inRow
            newRow = firstPart V.++ (value `V.cons` remainderPart)
        in 
        if L.null modList then (newRow, modList)
        else 
            let (nextRowNumber, nextColNumber, nextValue) = L.head modList
            in
            if nextRowNumber /= rowNumber then (newRow, modList) 
            else modifyRow newRow nextColNumber nextValue rowNumber (L.tail modList)
            




