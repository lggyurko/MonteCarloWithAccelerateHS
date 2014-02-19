% include lhs2TeX.sty
%include agda.fmt
%include lhs2TeX.fmt
% include polycode.fmt 



> {-# LANGUAGE UnicodeSyntax #-}
> module Main (main) where
> import MonteCarloExample as MC
> import Mrg32k3a
> import Data.Array.Accelerate as A
> import Data.Array.Accelerate.CUDA as I 
> import Unicode


>
> main :: IO()
> main =  do 
>    putStr "Average value and standard deviation of error = " 
>    putStr $ show mc1
>    putStr ", "
>    putStrLn $ show stdev 

> type Sc = Float -- Float
> generateNormalP = generateNormalPairF -- generateNormalPairF

Initialising the parameters. 

> n = 1000 :: Int
> m = 100 :: Int

Parameters of the BS model.

> r       = 0.05 :: Sc
> sigma   = 0.1 :: Sc
> t       = 1.0 :: Sc
> s0      = 1.0 :: Sc
> rho     = 0.5 :: Sc
> k       = 1.0 :: Sc

Some precomputed constants.

> v10 = (0.0,0.0) :: (Sc,Sc)
> rm = 1.0 / Prelude.fromIntegral m :: Sc
> erm = constant rm :: Exp Sc
> rn = 1.0 / Prelude.fromIntegral n :: Sc
> ern = constant rn :: Exp Sc
> rnm = rn * rm
> ernm = ern * erm
> rsq2m = sqrt $ 0.5 * rm
> ersq2m = constant rsq2m :: Exp Sc
> ek = constant k :: Exp Sc
> ekk = constant 0.1 :: Exp Sc

> er      = constant r :: Exp Sc
> esigma  = constant sigma :: Exp Sc
> erho    = constant rho :: Exp Sc
> et      = constant t :: Exp Sc
> edt     = et * erm :: Exp Sc
> esqdt   = sqrt edt  :: Exp Sc
> ealpha  = sqrt $ 1.0 - erho * erho :: Exp Sc
> econ1   = er * edt + 1.0 :: Exp Sc
> econ2   = esqdt * esigma :: Exp Sc 
> ezero   = constant 0 :: Exp Sc
> edisc   = exp $ er * et * (-1.0) :: Exp Sc

Actually, there is no need to take substeps when the underlying stock price follows a geometric Brownian motion and we only need the price at maturity. Here, we do that for testing purposes only. 

> gbm2dEulerMaruyama ::  Exp (Sc, Sc) → Exp (Sc,Sc) 
>                        → Exp (Sc,Sc)
> gbm2dEulerMaruyama sn zs = 
>    let
>       (s1,s2) = unlift sn :: (Exp Sc, Exp Sc)
>       (z1,z2) = unlift zs :: (Exp Sc, Exp Sc)
>       z2' = erho * z1 + ealpha * z2
>       s1' = s1 * (econ1 + econ2*z1)
>       s2' = s2 * (econ1 + econ2*z2')
>    in 
>       lift (s1',s2')

> payoff :: Exp (Sc, Sc) → Exp (Sc, Sc)
> payoff sn = 
>    let
>       (s1, s2) = unlift sn :: (Exp Sc, Exp Sc)
>       v =  ((abs (s1 - ek)) <* ekk &&* (abs (s2 - ek)) <* ekk) 
>            ? (edisc,ezero) :: Exp Sc
>    in 
>       lift (ern * v, ern * v * v) 

A particular aggregate function:

> sumP :: Exp (Sc,Sc) → Exp (Sc,Sc) → Exp (Sc,Sc)
> sumP p1 p2 =
>    let
>       (x1,y1) = unlift p1 :: (Exp Sc, Exp Sc)
>       (x2,y2) = unlift p2 :: (Exp Sc, Exp Sc)
>    in 
>       lift (x1 + x2, y1 + y2)

> mc = simpleMC n m state (s0, s0) generateNormalP gbm2dEulerMaruyama payoff sumP
>
> (mc1, mc2) = indexArray (I.run mc) Z
> stdev = sqrt $ rn * (mc2 - mc1 * mc1)




