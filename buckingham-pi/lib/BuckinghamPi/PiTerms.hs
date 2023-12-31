module BuckinghamPi.PiTerms where

-- import the following two libraries for the matrix data type, and matrix
-- functions
import qualified Data.Matrix as DM (Matrix, transpose, fromLists, toLists, detLU, zero, rref,
                                    (<|>), ncols, colVector, getCol, submatrix, (<->), scaleMatrix,
                                    identity, luDecomp', multStd, matrix, inverse)
import qualified Numeric.Matrix as NM (Matrix, fromList, toList, rank, zero, inv, isEmpty, isZero)
-- import Ratio library to handle the Rational data type
import Data.Ratio
-- import Maybe library to handle the Maybe data type
import Data.Maybe
-- import List library for extra capability in list permutations, as well as
-- importing two more list functions for list searching and testing uniqueness
import Data.List
import Data.List.Split (endBy)
import Data.List.Unique (isUnique)


-- create the Dim data type to represent dimensions of physical quantities
-- where there are a selection of 26 dimensions (A to Z)
data Dim = A | B | C | D | E | F | G | H | I | J | K | L | M | N | O | P | Q | R | S | T | U | V | W | X | Y | Z
  deriving (Eq, Show, Ord, Read)

-- define type synonym Pow to represent the power of a dimension 
type Pow = Rational
-- define type synonym Var to represent the name of a physical variable 
type Var = String
-- define type synonym DimQu to represent a dimensional quantity
type DimQu = [(Dim,Pow)]
-- define a type synonym to represent a dimensional matrix
type DimMat = DM.Matrix Pow
-- define a type synonym to represent a kernel matrix containing
-- the kernel vectors of a dimensional matrix
type Kernel = DM.Matrix Pow

-- the above two type synonyms imply that the DM matrix type is the default type


-- ###################################################
-- # Helpful "gluing" functions that are fundamental #
-- ###################################################
-- a function to sort a dimensional quantity, DimQu, so that dimensions
-- in [(dim,pow)] are not repeated, but added/subtracted together
sortDimQu' :: DimQu -> DimQu
sortDimQu' ((dim,pow):[]) = [(dim,pow)]
sortDimQu' ((dim1,pow1):(dim2,pow2):otherdims)
  | dim1 == dim2 = sortDimQu' ((dim1,pow1 + pow2):otherdims)
  | dim1 /= dim2 = (dim1,pow1) : sortDimQu' ((dim2,pow2):otherdims)
-- a function to sort a dimensional quantity, DimQu, including a preprocessing
-- of the input so all dimensions are in alphabetical order
sortDimQu :: DimQu -> DimQu
sortDimQu = sortDimQu' . sort

-- a function to sort a list of dimensional quantities, [DimQu]
sortDimQus :: [DimQu] -> [DimQu]
sortDimQus = map sortDimQu

-- a function that multiplies two dimensional quantities (meaning the powers of
-- their common dimensions are added together)
(|*|) :: DimQu -> DimQu -> DimQu
(|*|) dimqu_1 dimqu_2 = sortDimQu $ dimqu_1 ++ dimqu_2

-- a function that multiplies two dimensional quantities (meaning the powers of
-- their common dimensions are added together)
(|/|) :: DimQu -> DimQu -> DimQu
(|/|) dimqu_1 dimqu_2 = sortDimQu $ dimqu_1 ++ (map (\(dim,pow) -> (dim,(-pow))) dimqu_2)

-- function to create a column of a dimensional matrix as a list of powers of
-- all separate dimensions in a dimensional quantity
makeDimMatCol :: DimQu -> [Pow]
makeDimMatCol [] = []
makeDimMatCol ((dim,pow):otherdims) = pow : makeDimMatCol otherdims

-- function to make a dimensional matrix from a list of dimensional quantities
makeDimMat :: [DimQu] -> DimMat
makeDimMat = DM.transpose . DM.fromLists . map makeDimMatCol . sortDimQus

-- function to convert a matrix from the DM library to a matrix in the NM
-- library's data type
toNM :: DM.Matrix Pow -> NM.Matrix Pow
toNM = NM.fromList . DM.toLists

-- function to convert a matrix from the NM library to a matrix in the DM
-- library's data type
toDM :: NM.Matrix Pow -> DM.Matrix Pow
toDM = DM.fromLists . NM.toList

-- function to convert a RealFrac number to an integer by rounding it
toInt n = (round n) :: Int 

-- function to get the rank of a (dimensional) matrix
getRank :: DM.Matrix Pow -> Int
getRank =  toInt . NM.rank . toNM

-- function to get the determinant of a (dimensional) matrix
getDet :: DM.Matrix Pow -> Int
getDet = toInt . DM.detLU

-- function to get the inverse of a matrix
getInverse :: DM.Matrix Pow -> Maybe (DM.Matrix Pow)
getInverse = (either (\_ -> Nothing) Just) . DM.inverse

-- function to get the reduced row echelon form of a matrix
getRREF :: DM.Matrix Pow -> Maybe (DM.Matrix Pow)
getRREF = (either (\_ -> Nothing) Just) . DM.rref

-- function to find the unique combinations of elements from an original list
-- where the unique combinations are only of a certain length, n
-- example: combinations 2 "abc" == ["ab","ac","bc"]
combinations :: Int -> [a] -> [[a]]
combinations n = filter ((== n) . length) . subsequences

-- function to get the number of pi terms expected from a dimensional matrix
getNPiTerms :: DimMat -> Int
getNPiTerms dimMat = n - k
  where n = DM.ncols dimMat
        k = getRank dimMat

-- function to check if a matrix is composed of all zero elements
isMatZero :: DM.Matrix Pow -> Bool
isMatZero = NM.isZero . toNM


-- #####################
-- #  LU Decomposition #
-- #####################

-- LU decomposition is needed so that the number of pi terms of a dimensional
-- matrix can be obtained even when that dimensional matrix may have column
-- configurations where the whole matrix is rank deficient. This assumes that
-- getRank can only work on square matrix segments of rectangular matrices

-- function to get LU decomnposition with full pivoting
-- (where the output is (U,L,P,Q,d,e):
--     U is an upper triangular matrix.
--     L is an unit lower triangular matrix.
--     P,Q are permutation matrices.
--     d,e are the determinants of P and Q respectively.
--     PMQ = LU)
getLUDecomp :: DM.Matrix Pow -> (DM.Matrix Pow, DM.Matrix Pow, DM.Matrix Pow, DM.Matrix Pow, Pow, Pow)
getLUDecomp = fromMaybe (DM.zero 1 1,DM.zero 1 1,DM.zero 1 1,DM.zero 1 1,0,0) . DM.luDecomp'

-- function to only obtain the 
getLUDecompU :: DM.Matrix Pow -> DM.Matrix Pow
getLUDecompU = (\(u,_,_,_,_,_) -> u) . getLUDecomp


-- ####################################################
-- # Method of Repeating Variables (Helper Functions) #
-- ####################################################

-- function to get a list of unique combinations of non repeating variables
-- from a list
getNonRepVarsList :: Int -> [a] -> [[(Int,a)]]
getNonRepVarsList n ls = combinations n (zip [1..nLs] ls)
  where nLs = length ls

-- function to get a list of all unique combinations of non-repeating variables
getRepVarsCombsPM :: Eq a => Int -> [a] -> [[a]]
getRepVarsCombsPM n dimQus = map (snd . unzip) [(dimQusNumd \\ nonRepVars) ++ nonRepVars | nonRepVars <- nonRepVarsListNumd]
  where nonRepVarsListNumd = getNonRepVarsList n dimQus
        dimQusNumd = zip [1..(length dimQus)] dimQus

-- function to split up a dimensional matrix into a matrix of repeating
-- variables (called X, or P), and a matrix of non-repeating variables
-- (called Y, or Q)
getXY :: DM.Matrix Pow -> (DM.Matrix Pow,DM.Matrix Pow)
getXY dimMat = (x,y)
  where rank = getRank dimMat
        nCols = DM.ncols dimMat
        x = DM.submatrix 1 rank 1 rank $ dimMat
        y = DM.submatrix 1 rank (rank + 1) nCols $ dimMat


-- ####################################################
-- #  Method of Repeating Variables (Main Functions)  #
-- ####################################################

-- function to get all unique combinations of repeating variables and
-- non-repeating variables
getRepVarsCombs :: ([Var],[DimQu]) -> [([Var],[DimQu])]
getRepVarsCombs (vars,dimQus) = zip varsCombs dimQusCombs
  where nNonRepVars = getNPiTerms . getLUDecompU . makeDimMat $ dimQus
        dimQusCombs = getRepVarsCombsPM nNonRepVars dimQus :: [[DimQu]]
        varsCombs = getRepVarsCombsPM nNonRepVars vars :: [[Var]]

-- function to get the dimensional matrices from the list of all
-- unique combinations of repeating/non-repeating variables
getDimMatFromDimQus :: [([Var],[DimQu])] -> [([Var],DimMat)]
getDimMatFromDimQus varsDimQusList = (map . fmap) makeDimMat varsDimQusList

-- function to filter out dimensional matrices that are singular due to their
-- column configurations
getNonSingDimMats :: [([Var],DimMat)] -> [([Var],DimMat)]
getNonSingDimMats varsDimMatList = filter filterFunc varsDimMatList
  where filterFunc = \(_,dimMat) -> if (getDet dimMat /= 0) then True else False

-- function to filter out all dimensional matrices that would yield duplicate
-- pi terms via filtering out dimensional matrices that have the same RREF
getNonDupVarsDimQus :: [([Var],DimMat)] -> [([Var],DimMat)]
getNonDupVarsDimQus varsDimMatList = nonDupVarsRREFDimMatList
  where rrefDimMats = map getRREF . snd . unzip $ varsDimMatList
        varsDimMatRREFList = (map . fmap) ((fromMaybe (DM.zero 1 1)) . getRREF) varsDimMatList
        nonDupRREFDimMats = map (fromMaybe (DM.zero 1 1)) . filter isJust $ nub rrefDimMats
        nonDupVarsRREFDimMatList =
          [fromMaybe ([""],DM.zero 1 1) . find ((== nonDupRREFDimMat) . ((fromMaybe (DM.zero 1 1)) . getRREF) . snd) $ varsDimMatList | nonDupRREFDimMat <- nonDupRREFDimMats]

-- function to obtain the kernel matrix from all dimensional matrices using the
-- equation v = -inv(P)*Q*e or v = -inv(X)*Y*e where v and e constitute the
-- unknown and known (canonical basis vector matrix, e) parts of the kernel
-- matrix
getKernels :: [([Var],DimMat)] -> [([Var],Kernel)]
getKernels varsDimMatList = (map . fmap) getKernel varsDimMatList
  where getKernel dimMat = (DM.scaleMatrix (-1) (DM.multStd xInv (DM.multStd y canVecBasMat))) DM.<-> canVecBasMat
          where dimMatRREF = fromMaybe (DM.zero 1 1) $ getRREF dimMat
                (x,y) = getXY dimMatRREF
                xInv = fromMaybe (DM.zero 1 1) $ getInverse x
                canVecBasMat = DM.identity $ getNPiTerms dimMatRREF

-- a function wrap all main functions together in sequence
getSetsOfPiTerms :: ([Var],[DimQu]) -> [([Var],Kernel)]
getSetsOfPiTerms = getKernels .
                   getNonDupVarsDimQus .
                   getNonSingDimMats .
                   getDimMatFromDimQus .
                   getRepVarsCombs
