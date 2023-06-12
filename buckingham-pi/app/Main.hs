module Main where

import BuckinghamPi.TextParsing (getSetsOfPiTermsLaTeXStrings)
import System.Environment
import System.IO

main :: IO ()
main = do
  args <- getArgs
  case args of
    [inFile] -> do varsDimQusString <- readFile inFile
                   putStrLn $ getSetsOfPiTermsLaTeXStrings varsDimQusString
    _ -> putStrLn "ERROR: Wrong number of arguments (expected one)"
  
