% include lhs2TeX.sty
%include agda.fmt
%include lhs2TeX.fmt
% include polycode.fmt 

> {-# LANGUAGE UnicodeSyntax #-}
> module BlackScholes (callOption) where
> import Data.Number.Erf
> import Unicode


> rsq2 = sqrt 0.5 :: Double
>
> callOption :: Double → Double → Double → Double → Double → Double
> callOption s k r t sigma = 
>   let
>      sqv = sigma * (sqrt t)
>      v = sqv * sqv
>      d1 = ((log (s / k)) + r * t + 0.5 * v) / sqv
>      d2 = d1 - sqv
>   in
>      s * (normcdf d1) - k * (normcdf d2) * (exp $ (-1.0) *  r * t)