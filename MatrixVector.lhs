%include formats.fmt
%include agda.fmt
%include lhs2TeX.fmt


This module defines data types and operations to model and update the state of the random number generator, and to extract value from the state. 

> {-# LANGUAGE UnicodeSyntax #-}
> module MatrixVector 
> -- (prodEMMP, prodEMVP, N, Triple, State, Mat, MatPair) 
>   where
> import Unicode 
> import Data.Array.Accelerate as A


The following data types model vectors, states and 3-by-3 matrices respectively. 
The particular form ensures that these objects can be used as element types in Accelerate Arrays. 

> type N = Int64
> type Triple e   = (e,e,e)
> type State e    = (Triple e, Triple e)
> type Mat e      = Triple (Triple e)
> type MatPair e  = (Mat e, Mat e) 

The key operations of this module are the paired versions of the matrix-matrix multiplication and the matrix-vector multiplication on GPUs. 

Paired version of the matrix-matrix multiplication. 

> prodEMMP :: Exp N → Exp N → Exp (MatPair N) → Exp (MatPair N) → Exp (MatPair N)
> prodEMMP m1 m2 a b = 
>    let
>       (a1,a2) = unliftMP a :: MatPair (Exp N)
>       (b1,b2) = unliftMP b :: MatPair (Exp N)
>       uab1 = prodMM m1 a1 b1 :: Mat (Exp N)
>       uab2 = prodMM m2 a2 b2 :: Mat (Exp N)
>    in
>       liftMP (uab1, uab2)   


Paired version of the matrix-vector multiplication on GPUs:

> prodEMVP :: Exp N → Exp N → Exp (MatPair N) → Exp (State N) → Exp (State N)
> prodEMVP m1 m2 a b = 
>    let
>       (a1, a2) = unliftMP a :: MatPair (Exp N)
>       (b1, b2) = unliftS b :: State (Exp N)
>       uab1 = prodMV m1 a1 b1 :: Triple (Exp N)
>       uab2 = prodMV m2 a2 b2 :: Triple (Exp N)
>    in
>       liftS (uab1, uab2)

% ----------------------------------------------------
% ----------------------------------------------------
The rest of the module contains a few helper functions. 

> tr1 :: Triple e → e
> tr1 (a,_,_) = a
> tr2 :: Triple e → e
> tr2 (_,a,_) = a
> tr3 :: Triple e → e
> tr3 (_,_,a) = a

Modulo-product. We assume that the module is smaller than $2^{32}$, hence we work elements smaller than 2^32. 

> pr :: Exp N → Exp N → Exp N → Exp N
> pr m a b = 
>    let
>       mod' x = mod x m
>       a1 = div a 2
>       a2 = a - a1
>       b1 = div b 2
>       b2 = b - b1
>       ab11 = mod' $ a1 * b1
>       ab12 = mod' $ a1 * b2
>       ab21 = mod' $ a2 * b1
>       ab22 = mod' $ a2 * b2
>    in 
>       mod' $ ab11 + ab12 + ab21 + ab22

Dot-product:

> prodDot :: Exp N → Triple (Exp N) → Triple (Exp N) → Exp N
> prodDot m a b = 
>    let
>       mod' x = mod x m
>       pr' x y = pr m x y
>       (a1,a2,a3) = a 
>       (b1,b2,b3) = b 
>    in 
>       mod' $ (pr' a1 b1) + (pr' a2 b2) + (pr' a3 b3)

Matrix-vector multiplication:

> prodMV :: Exp N → Mat (Exp N) → Triple (Exp N) → Triple (Exp N)
> prodMV m a b = 
>    let
>       (a1,a2,a3) = a 
>       c1 = prodDot m a1 b
>       c2 = prodDot m a2 b
>       c3 = prodDot m a3 b
>    in 
>      (c1, c2, c3)   

Matrix-matrix multiplication:

> prodMM :: Exp N → Mat (Exp N) → Mat (Exp N) → Mat (Exp N)
> prodMM m a b = 
>    let
>       ((b11,b12,b13),(b21,b22,b23),(b31,b32,b33)) = b 
>       (c11,c21,c31) = prodMV m a (b11, b21, b31)                             
>       (c12,c22,c32) = prodMV m a (b12, b22, b32)
>       (c13,c23,c33) = prodMV m a (b13, b23, b33)     
>   in
>     ((c11,c12,c13),(c21,c22,c23),(c31,c32,c33))

Matrix-matrix multiplication on GPUs:

> prodEMM :: Exp N → Exp (Mat N) → Exp (Mat N) → Exp (Mat N)
> prodEMM m a b = 
>    let
>       ua = unliftM a
>       ub = unliftM b
>       uab = prodMM m ua ub
>    in
>       liftM uab   

Matrix-vector multiplication on GPUs:

> prodEMV :: Exp N → Exp (Mat N) → Exp (Triple N) → Exp (Triple N)
> prodEMV m a b = 
>    let
>       ua = unliftM a :: Mat (Exp N)
>       ub = unliftT b :: Triple (Exp N)
>       uab = prodMV m ua ub :: Triple (Exp N)
>    in
>       liftT uab   

Lift and unlift operations are required when manipulating vectors, matrices and states element-wise on the GPU. 

> liftT :: (Elt e) ⇒ Triple (Exp e) → Exp (Triple e)
> liftT = lift 

> liftS :: (Elt e) ⇒ State (Exp e) → Exp (State e)
> liftS a = 
>    let
>       (a1, a2) = a
>    in 
>       lift (liftT a1, liftT a2)

> liftM :: (Elt e) ⇒ Mat (Exp e) → Exp (Mat e)
> liftM a =
>    let
>       (a1, a2, a3) = a
>    in 
>       liftT (liftT a1, liftT a2, liftT a3)

> liftMP :: (Elt e) ⇒ MatPair (Exp e) → Exp (MatPair e)
> liftMP a =
>    let
>       (a1, a2) = a
>    in 
>       lift (liftM a1, liftM a2)

> unliftT :: (Elt e) ⇒ Exp (Triple e) → Triple (Exp e)
> unliftT = unlift

> unliftS :: (Elt e) ⇒ Exp (State e) → State (Exp e)
> unliftS a =   
>    let
>       (ea1, ea2) = unlift a 
>    in 
>       (unliftT ea1, unliftT ea2)

> unliftM :: (Elt e) ⇒ Exp (Mat e) → Mat (Exp e)
> unliftM a =   
>    let
>       (ea1, ea2, ea3) = unliftT a 
>    in 
>       (unliftT ea1, unliftT ea2, unliftT ea3)

> unliftMP :: (Elt e) ⇒ Exp (MatPair e) → MatPair (Exp e)
> unliftMP a =   
>    let
>       (ea1, ea2) = unlift a 
>    in 
>       (unliftM ea1, unliftM ea2)
