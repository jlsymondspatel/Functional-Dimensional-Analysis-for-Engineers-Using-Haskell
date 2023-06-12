module BuckinghamPi.PiTerms where

import qualified Data.Matrix as DM (Matrix, transpose, fromLists, toLists, detLU, zero, rref,
                                    (<|>), ncols, colVector, getCol, submatrix, (<->), scaleMatrix,
                                    identity, luDecomp', multStd, matrix, inverse)
import qualified Numeric.Matrix as NM (Matrix, fromList, toList, rank, zero, inv, isEmpty, isZero)
import Data.Ratio
import Data.Maybe
import Data.List
import Data.List.Split (endBy)
import Data.List.Unique (isUnique)


data Dim = A | B | C | D | E | F | G | H | I | J | K | L | M | N | O | P | Q | R | S | T | U | V | W | X | Y | Z
  deriving (Eq, Show, Ord, Read)


type Pow = Rational
type Var = String
type DimQu = [(Dim,Pow)]
type DimMat = DM.Matrix Pow
type Kernel = DM.Matrix Pow


-- ###################################################
-- # Helpful "gluing" functions that are fundamental #
-- ###################################################

sortDimQu' :: DimQu -> DimQu
sortDimQu' ((dim,pow):[]) = [(dim,pow)]
sortDimQu' ((dim1,pow1):(dim2,pow2):otherdims)
  | dim1 == dim2 = sortDimQu' ((dim1,pow1 + pow2):otherdims)
  | dim1 /= dim2 = (dim1,pow1) : sortDimQu' ((dim2,pow2):otherdims)

sortDimQu :: DimQu -> DimQu
sortDimQu = sortDimQu' . sort

sortDimQus :: [DimQu] -> [DimQu]
sortDimQus = map sortDimQu

(|*|) :: DimQu -> DimQu -> DimQu
(|*|) dimqu_1 dimqu_2 = sortDimQu $ dimqu_1 ++ dimqu_2

(|/|) :: DimQu -> DimQu -> DimQu
(|/|) dimqu_1 dimqu_2 = sortDimQu $ dimqu_1 ++ (map (\(dim,pow) -> (dim,(-pow))) dimqu_2)

makeDimMatCol :: DimQu -> [Pow]
makeDimMatCol [] = []
makeDimMatCol ((dim,pow):otherdims) = pow : makeDimMatCol otherdims

makeDimMat :: [DimQu] -> DM.Matrix Pow
makeDimMat = DM.transpose . DM.fromLists . map makeDimMatCol . sortDimQus

toNM :: DM.Matrix Pow -> NM.Matrix Pow
toNM = NM.fromList . DM.toLists

toDM :: NM.Matrix Pow -> DM.Matrix Pow
toDM = DM.fromLists . NM.toList

toInt n = (round n) :: Int 

getRank :: DM.Matrix Pow -> Int
getRank =  toInt . NM.rank . toNM

getDet :: DM.Matrix Pow -> Int
getDet = toInt . DM.detLU

getInverse :: DM.Matrix Pow -> Maybe (DM.Matrix Pow)
getInverse = (either (\_ -> Nothing) Just) . DM.inverse

getRREF :: DM.Matrix Pow -> Maybe (DM.Matrix Pow)
getRREF = (either (\_ -> Nothing) Just) . DM.rref

combinations :: Int -> [a] -> [[a]]
combinations n = filter ((== n) . length) . subsequences

getNPiTerms :: DM.Matrix Pow -> Int
getNPiTerms dimMat = n - k
  where n = DM.ncols dimMat
        k = getRank dimMat

isMatZero :: DM.Matrix Pow -> Bool
isMatZero = NM.isZero . toNM


-- #####################
-- #  LU Decomposition #
-- #####################

-- function to get LU decomnposition with full pivoting
-- (where the output is (U,L,P,Q,d,e):
--     U is an upper triangular matrix.
--     L is an unit lower triangular matrix.
--     P,Q are permutation matrices.
--     d,e are the determinants of P and Q respectively.
--     PMQ = LU)
getLUDecomp :: DM.Matrix Pow -> (DM.Matrix Pow, DM.Matrix Pow, DM.Matrix Pow, DM.Matrix Pow, Pow, Pow)
getLUDecomp = fromMaybe (DM.zero 1 1,DM.zero 1 1,DM.zero 1 1,DM.zero 1 1,0,0) . DM.luDecomp'

getLUDecompU :: DM.Matrix Pow -> DM.Matrix Pow
getLUDecompU = (\(u,_,_,_,_,_) -> u) . getLUDecomp


-- ####################################################
-- # Method of Repeating Variables (Helper Functions) #
-- ####################################################

getNonRepVarsList :: Int -> [a] -> [[(Int,a)]]
getNonRepVarsList n ls = combinations n (zip [1..nLs] ls)
  where nLs = length ls

getRepVarsCombsPM :: Eq a => Int -> [a] -> [[a]]
getRepVarsCombsPM n dimQus = map (snd . unzip) [(dimQusNumd \\ nonRepVars) ++ nonRepVars | nonRepVars <- nonRepVarsListNumd]
  where nonRepVarsListNumd = getNonRepVarsList n dimQus
        dimQusNumd = zip [1..(length dimQus)] dimQus

getXY :: DM.Matrix Pow -> (DM.Matrix Pow,DM.Matrix Pow)
getXY dimMat = (x,y)
  where rank = getRank dimMat
        nCols = DM.ncols dimMat
        x = DM.submatrix 1 rank 1 rank $ dimMat
        y = DM.submatrix 1 rank (rank + 1) nCols $ dimMat


-- ####################################################
-- #  Method of Repeating Variables (Main Functions)  #
-- ####################################################

getRepVarsCombs :: ([Var],[DimQu]) -> [([Var],[DimQu])]
getRepVarsCombs (vars,dimQus) = zip varsCombs dimQusCombs
  where nNonRepVars = getNPiTerms . getLUDecompU . makeDimMat $ dimQus
        dimQusCombs = getRepVarsCombsPM nNonRepVars dimQus :: [[DimQu]]
        varsCombs = getRepVarsCombsPM nNonRepVars vars :: [[Var]]

getDimMatFromDimQus :: [([Var],[DimQu])] -> [([Var],DimMat)]
getDimMatFromDimQus varsDimQusList = (map . fmap) makeDimMat varsDimQusList

getNonSingDimMats :: [([Var],DimMat)] -> [([Var],DimMat)]
getNonSingDimMats varsDimMatList = filter filterFunc varsDimMatList
  where filterFunc = \(_,dimMat) -> if (getDet dimMat /= 0) then True else False

getNonDupVarsDimQus :: [([Var],DimMat)] -> [([Var],DimMat)]
getNonDupVarsDimQus varsDimMatList = nonDupVarsRREFDimMatList
  where rrefDimMats = map getRREF . snd . unzip $ varsDimMatList
        varsDimMatRREFList = (map . fmap) ((fromMaybe (DM.zero 1 1)) . getRREF) varsDimMatList
        nonDupRREFDimMats = map (fromMaybe (DM.zero 1 1)) . filter isJust $ nub rrefDimMats
        nonDupVarsRREFDimMatList =
          [fromMaybe ([""],DM.zero 1 1) . find ((== nonDupRREFDimMat) . ((fromMaybe (DM.zero 1 1)) . getRREF) . snd) $ varsDimMatList | nonDupRREFDimMat <- nonDupRREFDimMats]

getKernels :: [([Var],DimMat)] -> [([Var],Kernel)]
getKernels varsDimMatList = (map . fmap) getKernel varsDimMatList
  where getKernel dimMat = (DM.scaleMatrix (-1) (DM.multStd xInv (DM.multStd y canVecBasMat))) DM.<-> canVecBasMat
          where dimMatRREF = fromMaybe (DM.zero 1 1) $ getRREF dimMat
                (x,y) = getXY dimMatRREF
                xInv = fromMaybe (DM.zero 1 1) $ getInverse x
                canVecBasMat = DM.identity $ getNPiTerms dimMatRREF

getSetsOfPiTerms :: ([Var],[DimQu]) -> [([Var],Kernel)]
getSetsOfPiTerms = getKernels .
                   getNonDupVarsDimQus .
                   getNonSingDimMats .
                   getDimMatFromDimQus .
                   getRepVarsCombs
