module BuckinghamPiProbProp.ProbProp where

import qualified Statistics.Distribution as SD
import qualified Statistics.Distribution.Normal as SDN
import qualified Statistics.Distribution.Uniform as SDU
import Control.Applicative


data YS a = YS a
  deriving (Show,Ord,Eq)

instance Functor YS where
  fmap func (YS a) = YS (func a)

instance Applicative YS where
  pure = YS
  YS func <*> YS x = YS (func x)

data CDF a = CDF a
  deriving (Show,Ord,Eq)

instance Functor CDF where
  fmap func (CDF a) = CDF (func a)

instance Applicative CDF where
  pure = CDF
  CDF func <*> CDF x = CDF (func x)

data PDF a = PDF a
  deriving (Show,Ord,Eq)

instance Functor PDF where
  fmap func (PDF a) = PDF (func a)

instance Applicative PDF where
  pure = PDF
  PDF func <*> PDF x = PDF (func x)

fsCDFtoPDF (CDF a) = PDF a
fsPDFtoCDF (PDF a) = CDF a


type Transform = [Double] -> YS Double
type OriginalPDF = [Double] -> PDF Double

normpdf :: Double -> Double -> Double -> Double
normpdf mu sigma x = SD.density (SDN.normalDistr mu sigma) x

normstdpdf :: Double -> Double
normstdpdf =  SD.density SDN.standard

unifpdf :: Double -> Double -> Double -> Double
unifpdf a b x = SD.density (SDU.uniformDistr a b) x

dx1 = 0.05 :: Double
dx2 = 0.05 :: Double
dxs = [dx1,dx2]

x1s = [0,dx1..20]
x2s = [0,dx2..5]
xss = [x1s,x2s]

g :: Transform
g [x1,x2] = pure (x1 * x2)
g _ = pure 0

fPDF :: OriginalPDF
fPDF [x1,x2] = pure ((normpdf 10 1 x1) * (unifpdf 1 3 x2))
fPDF _ = pure 0

dy = YS 0.1

pushPDF :: OriginalPDF -> Transform -> YS Double -> YS (PDF Double)
pushPDF f g y = fmap fsCDFtoPDF $ (fmap . fmap) (/ ((\(YS u) -> u) dy)) dFy
  where pushCDF yval = pure $ fsPDFtoCDF $ fmap (product dxs *) . fmap sum . sequenceA $ fCartProdFilt :: YS (CDF Double)
          where fCartProdFilt = map f . filter (\x1x2s -> g x1x2s <= yval) $ sequence xss
        dFy = (liftA2 . liftA2) (-) (pushCDF (liftA2 (+) y dy)) (pushCDF y) :: YS (CDF Double)
