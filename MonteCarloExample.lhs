%include formats.fmt
%include agda.fmt
%include lhs2TeX.fmt



> {-# LANGUAGE UnicodeSyntax #-}
> {-# LANGUAGE ScopedTypeVariables #-}
> module MonteCarloExample (simpleMC, state) where
> import Mrg32k3a 
> import MatrixVector
> import Unicode 
> import Data.Array.Accelerate as A


The simpleMC function consists of the following steps:

\begin{itemise}

\item initialise the generator states and model states,

\item define the update actions and pipes them together,

\item aggregates the final model values across threads. 

\end{itemise}

The function takes the following input arguments:

\begin{itemise}

\item number of threads, 

\item half the number of steps per thread

\item initial generator state

\item initial model state

\item an (Exp Double → Exp Double → Exp e → Exp e) update function

\item an (Exp e → Exp a) payoff function

\item an (Exp a → Exp a → Exp a) aggregator function

\end{itemise}

Types: \cd{a} : random variable(s), \cd{b} model value, \cd{c} payoff value. 

> simpleMC ::  ∀ a b c . (Elt a, Elt b, Elt c) ⇒ Int → Int → State N → b
>              → (Exp (State N) → Exp (State N, a)) -- rnd var. gen fn
>              → (Exp b → Exp a → Exp b)  -- model update fn     
>              → (Exp b → Exp c) -- payoff fn
>              → (Exp c → Exp c → Exp c) -- aggregator fn
>              → Acc (Array DIM0 c)
> simpleMC n m s v gf mf pf af =
>    let
>       en = constant n :: Exp Int
>       em = constant m :: Exp Int
>       ev = constant v :: Exp b
>       es = constant s :: Exp (State N)
>
>       -- initial generator states, one for each thread
>       genSts :: Acc (Array DIM1 (State N)) 
>       genSts  = initStates (unit es) (unit en) (unit em)
>
>       -- initial model states, one for each thread
>       modSts :: Acc (Array DIM1 b) 
>       modSts  = fill (index1 en) ev 
>
>       -- zip generator states and model states
>       initSts :: Acc (Array DIM1 (State N,b)) 
>       initSts = A.zip genSts modSts
>
>       -- update function that links the rv generator fn to model update fn
>       updFn :: Exp (State N, b) → Exp (State N, b)
>       updFn   = updateFn gf mf 
>
>       -- mapping the update function through the state vector
>       step :: Acc (Array DIM1 (State N,b)) → Acc (Array DIM1 (State N,b))
>       step    = A.map updFn 
>
>       -- list of operations (do we need a list?)
>       steps :: [Acc (Array DIM1 (State N,b)) → Acc (Array DIM1 (State N,b))]
>       steps   = Prelude.replicate m step
>
>       -- step operations piped together
>       finalSts = Prelude.foldl1 (>->) steps initSts :: Acc (Array DIM1 (State N,b))
>
>       -- extract final model values
>       (_,finalVals) = A.unzip finalSts :: Acc (Array DIM1 b)
>
>       -- mapping payoff over final states
>       payoffVals = A.map pf finalVals :: Acc (Array DIM1 c)
>
>    in
>       -- aggregation of payoff values
>       fold1 af payoffVals :: Acc (Aray DIM0 c)


Initial state of the generator.

> s1 = (1062452522,2961816100,342112271) :: Triple N
> s2 = (2854655037,3321940838,3542344109) :: Triple N
> state = (s1, s2)

Initial state is taken from \url{http://www.math.purdue.edu/~lucier/srfi-27/mrg32k3a.scm}

Following \url{http://www.cse.unsw.edu.au/~tmcdonell/presentations/2013-lambdajam-workshop.pdf}, arguments that may change should be passed to kernel function in as Arrays. 
