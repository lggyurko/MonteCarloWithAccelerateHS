%include formats.fmt
%include agda.fmt
%include lhs2TeX.fmt


> {-# LANGUAGE ScopedTypeVariables #-}
> {-# LANGUAGE UnicodeSyntax #-}
>
> module Mrg32k3a 
>    ( initStates, generateNormalPair, generateNormalPairF, 
>      generateUniform, generateUniformF, updateFn,
>      stateUpdate, stateToUniform )  where
> import MatrixVector
> import Unicode 
> import Data.Array.Accelerate as A

The following function initialises $n$ states, from each one $m$ steps can be taken, then the stats will overlap. The state in the first thread is $s$, which determines all the other states. 

> initStates ::  Acc (Scalar (State N)) → Acc (Scalar Int) → Acc (Scalar Int) → 
>                Acc (Array DIM1 (State N))
> initStates s n m = initStatesGen s n m am1 am2 ac

In a typical Monte-Carlo application, we generate a bunch of initial generator state and model state pairs (one for each thread), then repeatedly map an update function over them. Finally, the final model states get aggregated. The update function is constructed from a random variable generator (takes generator state and returns the pair of the new state and a random variable generated from the input state), and a model state update function (takes a model state a random variable and returns the new state of the model state).

> updateFn ::  ∀ e a . (Elt e, Elt a) ⇒ 
>               (Exp (State N) → Exp (State N, a)) → -- gen state, new gen state, rnd var
>               (Exp e → Exp a → Exp e) → -- model state, rnd var, new model state
>                Exp (State N, e) → Exp (State N, e) -- gen state, mod state
> updateFn gf mf sv =
>    let
>       (s,v)  = unlift sv :: (Exp (State N), Exp e)
>       (s1,r) = unlift $ gf s :: (Exp (State N), Exp a)
>    in
>       lift (s1, mf v r) :: Exp (State N, e)

In the above function, \cd{f} models the function that takes two generator states and a model state, then returns the new model state. This is useful for example when we use the Box-Muller method to convert a pair of uniform random variables to a pair of standard normals. 

State update generates a new generator state from the input state. 

> stateUpdate :: Exp (State N) → Exp (State N)
> stateUpdate = prodEMVP m1 m2 ec

Particular distributions. 

From generator state to uniform(0,1):

> stateToUniform :: Exp (State N) → Exp Double
> stateToUniform s =
>    let
>       (s1,s2) = unliftS s
>       x1 = tr1 s1
>       x2 = tr1 s2
>       x  = A.fromIntegral $ mod (x1 + x2) m1 :: Exp Double 
>    in
>       x * rm1 

Float version of stateToUniform - matters on GPUs (?)

> stateToUniformF :: Exp (State N) → Exp Float
> stateToUniformF s =
>    let
>       (s1,s2) = unliftS s
>       x1 = tr1 s1
>       x2 = tr1 s2
>       x  = A.fromIntegral $ mod (x1 + x2) m1 :: Exp Float 
>    in
>       x * rm1f 

Generator: returns the new state of the generator and a uniform(0,1) variable. 

> generateUniform :: Exp (State N) → Exp (State N, Double)
> generateUniform s = let s1 = stateUpdate s in lift (s1, stateToUniform s)

Float version of generateUniform - matters on GPUs (?)

> generateUniformF :: Exp (State N) → Exp (State N, Float)
> generateUniformF s = let s1 = stateUpdate s in lift (s1, stateToUniformF s)


Box-Muller method:

> tpi = constant 2.0 * pi :: Exp Double
> tpif = constant 2.0 * pi :: Exp Float

> boxMuller :: Exp (State N) → Exp (State N) → Exp (Double, Double)
> boxMuller s1 s2 = 
>    let
>       y1 = stateToUniform s1
>       y2 = stateToUniform s2
>       l = sqrt $ (-2.0) * (log y1) :: Exp Double
>       p = tpi * y2
>       c = cos p
>       s = sin p
>    in
>      lift (l * c, l * s)

Float version of boxMuller:

> boxMullerF :: Exp (State N) → Exp (State N) → Exp (Float,Float)
> boxMullerF s1 s2 = 
>    let
>       y1 = stateToUniformF s1
>       y2 = stateToUniformF s2
>       l = sqrt $ (-2.0) * (log y1) :: Exp Float
>       p = tpif * y2
>       c = cos p
>       s = sin p
>    in
>      lift (l * c, l * s)


Generator: returns the new state of the generator (after two updates!), and a pair of standard normals. 

> generateNormalPair :: Exp (State N) → Exp (State N, (Double, Double)) 
> generateNormalPair s =
>    let
>       s1 = stateUpdate s
>       s2 = stateUpdate s1
>       ns = boxMuller s s1 
>    in 
>       lift (s2, ns)


Float version of generateNormalPair

> generateNormalPairF :: Exp (State N) → Exp (State N, (Float,Float)) 
> generateNormalPairF s =
>    let
>       s1 = stateUpdate s
>       s2 = stateUpdate s1
>       ns = boxMullerF s s1 
>    in 
>       lift (s2, ns)


Marsaglia's polar method is a bit more tricky due to its acceptance-rejection nature. It might be possible to implement it by adding a counter to the model state. The update function increments the counter when the pair of uniforms are accepted, otherwise it does not change the state. Moreover, on the threads where the counter has reached the target, the model states are not modified any more. Difficulty: what if on some threads we run out of generator states before meeting the model state target due to many rejections? 
Possible solution: keep the model states (counters inclusive), regenerate the generator states by skipping ahead as needed, like in the Mandelbrot set example:
\url{http://community.haskell.org/~simonmar/slides/cadarache2012/7%20-%20accelerate.pdf}. 



The above approach generates random numbers on the GPU and consumes them shortly after generation. To consider: implement a version that only generates random numbers (from various distributions) to an array on the GPU, that can be either copied to the host memory or further processed on the GPU. 

Some constants:

> a1 = 1403580 :: N
> b1 = 810728 :: N
> m1 = constant $ 2^32 - 209 :: Exp N
> am1 = unit m1 :: Acc (Scalar N)
> rm1 = 1.0 / A.fromIntegral m1 :: Exp Double
> rm1f = 1.0 / A.fromIntegral m1 :: Exp Float
> a2 = 527612 :: N
> b2 = 1370589 :: N
> m2 = constant $ 2^32 - 22853 :: Exp N
> am2 = unit m2 :: Acc (Scalar N)
> c  = (c1,c2) :: MatPair N
> ec  = constant (c1,c2) :: Exp (MatPair N)
> ac = unit ec :: Acc (Scalar (MatPair N))

> c1 = ((0,a1,-b1),(1,0,0),(0,1,0)) :: Mat N
> c2 = ((a2,0,-b2),(1,0,0),(0,1,0)) :: Mat N
> idm = ((1,0,0),(0,1,0),(0,0,1)) :: Mat N
> eidm = constant (idm,idm) :: Exp (MatPair N)

Initialise states for $n$ threads and $m$ steps per thread. 

> initStatesGen ::  Acc (Scalar (State N)) → Acc (Scalar Int) → Acc (Scalar Int) → 
>                   Acc (Scalar N) → Acc (Scalar N) → Acc (Scalar (MatPair N)) →
>                   Acc (Array DIM1 (State N))
> initStatesGen s n m m1 m2 c =
>    let 
>       arr1   = fill (index1 $ the m) $ the c :: Acc (Array DIM1 (MatPair N)) 
>       fn :: Exp (MatPair N) → Exp (MatPair N) → Exp (MatPair N)
>       fn a b = prodEMMP (the m1) (the m2) a b
>       am     = the $ fold fn eidm arr1   -- A^m pair as Exp (MatPair Int)
>                                          -- is fold parallel-optimised?
>       arr2   = fill (index1 $ (the n) - 1) am :: Acc (Array DIM1 (MatPair N)) 
>       arr3   = A.scanl fn eidm arr2 :: Acc (Array DIM1 (MatPair N)) 
>    in
>       A.map (\ a → prodEMVP (the m1) (the m2) a (the s)) arr3   

