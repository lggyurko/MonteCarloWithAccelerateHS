
%include agda.fmt
%include lhs2TeX.fmt

This piece of code enables the use of unicode characters. 

\begin{code}
{-# LANGUAGE UnicodeSyntax #-}

module Unicode
where
\end{code}

\begin{code}
infixr 2 ∨
(∨)  ∷  Bool → Bool → Bool
a ∨ b  =  a || b
\end{code}

\begin{code}
infixr 3 ∧
(∧)  ∷  Bool → Bool → Bool
a ∧ b  =  a && b
\end{code}

\begin{code}               
infix 4 ≤, ≥
(≤), (≥)  ∷  (Ord a) ⇒ a → a → Bool
a ≤ b  =  a <= b
a ≥ b  =  a >= b
\end{code}

\begin{code}
infixr 9 ·
(·)  ∷  (b → c ) → (a → b) → (a → c)
f · g  =  \ x → f (g x)
\end{code}

∷ ⇒ ∀ → ← ⋯ ∨ ∧ ≤ ≥ ·