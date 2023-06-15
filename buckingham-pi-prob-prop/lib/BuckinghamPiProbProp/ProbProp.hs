module BuckinghamPiProbProp.ProbProp where

-- import Statistics library in order to use PDF functions for several
-- probability distributions
import qualified Statistics.Distribution as SD
import qualified Statistics.Distribution.Normal as SDN
import qualified Statistics.Distribution.Uniform as SDU
-- import Applicative library in order to use liftA2 and the Applicative
-- data type
import Control.Applicative

-- define data type to represent the context of a "Y space"
data YS a = YS a
  deriving (Show,Ord,Eq)

instance Functor YS where
  fmap func (YS a) = YS (func a)

instance Applicative YS where
  pure = YS
  YS func <*> YS x = YS (func x)

-- define data type to represent the context of a "CDF value space"
data CDF a = CDF a
  deriving (Show,Ord,Eq)

instance Functor CDF where
  fmap func (CDF a) = CDF (func a)

instance Applicative CDF where
  pure = CDF
  CDF func <*> CDF x = CDF (func x)

-- define data type to represent the context of a "PDF value space"
data PDF a = PDF a
  deriving (Show,Ord,Eq)

instance Functor PDF where
  fmap func (PDF a) = PDF (func a)

instance Applicative PDF where
  pure = PDF
  PDF func <*> PDF x = PDF (func x)

-- utility functions to convert from and to a PDF/CDF value space
fsCDFtoPDF (CDF a) = PDF a
fsPDFtoCDF (PDF a) = CDF a

-- define type synonyms that define the funcitonal type for a transformation of
-- random variables, g
type Transform = [Double] -> YS Double
-- define type synonyms that define the funcitonal type for a known original PDF
type OriginalPDF = [Double] -> PDF Double

-- define known PDF functions of known random variables
normpdf :: Double -> Double -> Double -> Double
normpdf mu sigma x = SD.density (SDN.normalDistr mu sigma) x

normstdpdf :: Double -> Double
normstdpdf =  SD.density SDN.standard

unifpdf :: Double -> Double -> Double -> Double
unifpdf a b x = SD.density (SDU.uniformDistr a b) x

-- define the spacing of the domain of the known (joint) PDF
dx1 = 0.1 :: Double
dx2 = 0.1 :: Double
dx3 = 0.1 :: Double
dxs = [dx1,dx2,dx3]

-- define the ranges of each argument to the domain of the known (joint) PDF
x1s = [12,(12+dx1)..18] -- l
x2s = [7,(7+dx2)..13] -- g
x3s = [2,(2+dx3)..8] -- tau
xss = [x1s,x2s,x3s]

-- define the function that transform the known random variables
g :: Transform
g [x1,x2,x3] = pure ((x1**(-0.5)) * (x2**0.5) * x3)
g _ = pure 0

-- define the whole joint PDF of all known random variables which may or may
-- not be independent
fPDF :: OriginalPDF
fPDF [x1,x2,x3] = pure ((normpdf 15 1 x1) * (normpdf 9.8 1 x2) * (normpdf 5.4 1 x3))
fPDF _ = pure 0

dy = YS 0.1

-- define the pushforward distribution function that obtains a PDF value
-- of the transformed random variable Y
pushPDF :: OriginalPDF -> Transform -> YS Double -> YS (PDF Double)
pushPDF f g y = fmap fsCDFtoPDF $ (fmap . fmap) (/ ((\(YS u) -> u) dy)) dFy
  where pushCDF yval = pure $ fsPDFtoCDF $ fmap (product dxs *) . fmap sum . sequenceA $ fCartProdFilt :: YS (CDF Double)
          where fCartProdFilt = map f . filter (\xs -> g xs <= yval) $ sequence xss
        dFy = (liftA2 . liftA2) (-) (pushCDF (liftA2 (+) y dy)) (pushCDF y) :: YS (CDF Double)
