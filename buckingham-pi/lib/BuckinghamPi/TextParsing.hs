module BuckinghamPi.TextParsing where

import qualified Data.Matrix as DM (Matrix, transpose, toLists)
import Data.Ratio
import Data.Maybe
import Data.List
import Data.List.Split (endBy)
import Data.List.Unique (isUnique)

import BuckinghamPi.PiTerms

-- #########################
-- # Textual LaTeX Parsing #
-- #########################

getExplicitUnitPows :: String -> String
getExplicitUnitPows [c] = if elem c ['A'..'Z']
                          then c : '^' : '{' : '1' : '}' : []
                          else [c]
getExplicitUnitPows (x1:x2:xs) = if elem x1 ['A'..'Z'] && elem x2 ['A'..'Z']
                                 then x1 : '^' : '{' : '1' : '}' : getExplicitUnitPows (x2:xs)
                                 else x1 : getExplicitUnitPows (x2:xs)

getVarsAndDimsStrings :: [String] -> ([Var],[String])
getVarsAndDimsStrings strings
  | not . and $ map (elem '=') strings = (take nVars $ repeat "invalid",take nVars $ repeat "#")
  | otherwise = (stringsLHSs,stringsRHSs)
  where stringsLHSs = map (takeWhile (/= '=')) strings
        stringsRHSs = map (getExplicitUnitPows . tail . dropWhile (/= '=') . filter (/= ' ')) strings
        nVars = length strings

getAllDims :: [String] -> [Dim]
getAllDims = map (read . (: [])) . sort . nub . filter ((flip elem) ['A'..'Z']) . concat

getDimQu :: [Dim] -> String -> DimQu
getDimQu allDims stringDims = dimQuWithAllDims
  where varDims = map (read . (: [])) . filter ((flip elem) ['A'..'Z']) $ stringDims :: [Dim]
        varPowsStrings = map (filter (/= '{')) . endBy "}" . filter ((flip elem) (['0'..'9']++['{','}','-'])) $ stringDims
        varPows = map (\x -> read (x ++ " % 1") :: Rational) varPowsStrings
        varDimsPows = sortDimQu $ zip varDims varPows
        dimQuWithAllDims = zip allDims (map (\dim -> fromMaybe 0 $ lookup dim varDimsPows) allDims)

getDimQus :: [String] -> [DimQu]
getDimQus stringDims = dimQus
  where allDims = getAllDims stringDims
        dimQus = map (getDimQu allDims) stringDims

getCheckedDimQus :: [String] -> [DimQu]
getCheckedDimQus stringDims
  | and $ map (checkIsInputGood) stringDims = getDimQus stringDims
  | otherwise = []
  where checkIsInputGood stringDim = if (nDims == nPows &&
                                         nDims == nCarrots &&
                                         nDims * 2 == nParens &&
                                         nPows == nCarrots &&
                                         nPows * 2 == nParens &&
                                         nOddChars == 0)
                                     then True
                                     else False
          where nDims = length $ filter ((flip elem) ['A'..'Z']) stringDim
                nCarrots = length $ filter ((flip elem) ['^']) stringDim
                nParens = length $ filter ((flip elem) ['{','}']) stringDim
                nPows = length . filter (not . null) . map (filter (/= '{')) . endBy "}" . filter ((flip elem) (['0'..'9']++['{','}','-'])) $ stringDim
                nOddChars = length $ filter (not . (flip elem) (['A'..'Z']++['0'..'9']++['{','}','^','-'])) stringDim

getVarsAndDimQus :: [String] -> ([Var],[DimQu])
getVarsAndDimQus inputString = (vars,dimQus)
  where (vars,rhs) = getVarsAndDimsStrings inputString
        dimQus = getCheckedDimQus rhs

getVarsMatAsVarsPowss :: ([Var],DM.Matrix Pow) -> ([Var],[[Pow]])
getVarsMatAsVarsPowss (vars,mat) = (vars,powss)
  where powss = (DM.toLists . DM.transpose) mat

printPow :: Pow -> String
printPow n = if (denominator n == 1)
             then show (numerator n)
             else map repl nString
  where nString = show n
        repl n
          | n == '%' = '/'
          | otherwise = n

getPiTermLaTeX :: Int -> [Var] -> [Pow] -> String
getPiTermLaTeX n vars pows = "\\Pi_{" ++ (show n) ++ "} = " ++ varsPows
  where varsPowsList = transpose $ [vars,powStringsLaTeX]
        powStrings = map printPow pows
        powStringsLaTeX = map (\p -> "^{" ++ p ++ "}") powStrings
        filterFuncs = filter (not . isSuffixOf ["^{0}"])
                       . map (\[var,pow] -> if pow == "^{1}" then [var,""] else [var,pow])
                       . map (\[var,pow] -> if take 6 pow == "^{(-1)" then [var,("^{-1")++(drop 6 pow)] else [var,pow])
        varsPows = concat . concat . filterFuncs $ varsPowsList

getSetOfPiTermsLaTeX :: ([Var],Kernel) -> String
getSetOfPiTermsLaTeX varsMat = result ++ "\\\\"
  where (vars,powss) = getVarsMatAsVarsPowss varsMat
        nPiTerms = length powss
        nthPiTerms = [1..nPiTerms]
        result = intercalate "\\qquad" $ zipWith3 getPiTermLaTeX nthPiTerms (replicate nPiTerms vars) powss

getSetsOfPiTermsLaTeX :: [([Var],Kernel)] -> [String]
getSetsOfPiTermsLaTeX = map getSetOfPiTermsLaTeX

getSetsOfPiTermsLaTeXStrings :: String -> String
getSetsOfPiTermsLaTeXStrings varsDimQusString
  | null dimQus == True = "ERROR #1: Invalid input stream"
  | det == 0 || isMatZero u = "ERROR #2: Dimensional Matrix is singular."
  | nPiTerms <= 0 = "ERROR #3: There are no extractable Pi terms."
  | otherwise = unlines . getSetsOfPiTermsLaTeX . getSetsOfPiTerms . getVarsAndDimQus . lines $ varsDimQusString
  where (vars,dimQus) = getVarsAndDimQus . lines $ varsDimQusString
        dimMat = makeDimMat . sortDimQus $ dimQus
        det = getDet dimMat
        u = getLUDecompU dimMat
        nPiTerms = getNPiTerms u
