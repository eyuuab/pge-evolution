! (register-module! ./utilities)
! (import! &self utilities:general-helpers)

;; bind python's round function to rd 
! (bind! rd (py-atom round))

;; createGenotype of a specified length in the range[0,1] and rounded off 
;; Parameters:- 
;;          $len length of the genotype tuple
;;          $prec signinfican digits for round function
;; Example :- 
;;          (createGenotype 6 2)
;;          -- 1st run -- (0.82 0.78 0.91 0.61 0.92 0.5)
;;          -- 2nd run -- (0.79 0.58 0.19 0.01 0.03 0.18)
;;          -- 3rd run -- (0.1 0.31 0.04 0.86 0.74 0.18)

(: (createGenotype (-> Number Number Expression)))
(= (createGenotype $len $prec)
    (if (== $len 0)
        ()
        (let* 
            (   
                ($first (rd ((py-atom random.random )) $prec))
                ($c (createGenotype (- $len 1) $prec))
            )
            (cons-atom $first $c)
            )))

;; Proabilistic Context Free Grammar (PCFG) representaion
;;          &grammar space created
;; Grammar rules added to the space as (non-terminal ((derivation_option) $ordering_index))
;; Example:- 
;;          (start ((expr) 0)) -- the first, and in this case the last, derivation rule for the "start" non-terminal 
;;          (op (("*") 1)) -- the second derivation rule for the "op" non-terminal
;;          probability values will be appended to this representation as follows
;;          (start ((expr prob_val) 0))
;;          (op (("*" prob_val) 1))
;;          prob_val being the probablity assigned to the particular derivation rule

! (bind! &grammar (new-space))

! (add-atom &grammar (start ((expr) 0)))
! (add-atom &grammar (expr (((expr op expr)) 0)))
! (add-atom &grammar (expr ((var) 1)))
! (add-atom &grammar (op (("+") 0)))
! (add-atom &grammar (op (("*") 1)))
! (add-atom &grammar (op (("-") 2)))
! (add-atom &grammar (var ((x) 0)))
! (add-atom &grammar (var ((1.0) 1)))

;; getNt -- gets the non-terminals from the grammar space

(= (getNt $grammarSpace)
    (collapse (unique (match &grammar ($nt $rule) $nt)))
)
;; updateRule -- updates the grammar space with rules conatining probablity values and deletes the old representation
;;          

(= (updateRule $grammarSpace $nt $len)
    (let () (match $grammarSpace ($nt (($x) $n)) (add-atom $grammarSpace ($nt (($x (rd (/ 1 $len) 2)) $n))))
    
        (match $grammarSpace ($nt (($x) $n)) (remove-atom $grammarSpace ($nt (($x) $n))))
    )
)
;; assignInitialProbs applies the updateRule function for all the non-terminals 

(= (assignInitialProbs $grammarSpace $ntSet)
        (if (== $ntSet ())
            (empty)
            (let* 
                (   
                    (($nt $tail) (decons $ntSet))
                    ($len (sum (collapse (match $grammarSpace ($nt $x) 1))))
                    (() (updateRule $grammarSpace $nt $len))
                )
                (assignInitialProbs $grammarSpace $tail)
            ))) 

;; ************ initialize initial probabality for grammar derivation options ***********

(= (initP $grammarSpace)
    (let $ntSet (getNt $grammarSpace)
        (assignInitialProbs $grammarSpace $ntSet)
    )
)
! (initP &grammar)

;; *********************************************************************************************
;; *********************************************************************************************
;; 
;; getOption gets the correct expansion option based on cummulative probablity values
;; Parameteres:- 
;;               $genotypeProb -- probabilty of derviation options from the genotype list of probablity values
;;               $initialCumProb -- the initial cummulative probablity that starts as 0 and keeps being updated
;;               $index -- index of the expansion options for a given non-terminal for ordered exploration,  it starts as 0
;; Examples:- 
;;               (getOption &grammar start 0.61 0 0) -- expr             
;;               (getOption &grammar expr 0.5 0 0) --  (expr op expr)             
;;               (getOption &grammar expr 0.7 0 0) -- var             
            
(= (getOption $grammarSpace $nt $genotypeProb $initialCumProb $index)
    (match $grammarSpace ($nt (($x $prob) $index)) 
            (if (<= $genotypeProb (+ $prob $initialCumProb)) 
                $x 
                (let* 
                    (
                        ($cumProb (+ $prob $initialCumProb))
                        ($next (+ $index 1))
                    )
                    (getOption $grammarSpace $nt $genotypeProb $cumProb $next)
                ))))
;; Original function used to get the expansion options of a literal non-terminal (not an expression like (expr or expr))

(= (ntByt $grammarSpace $symbol $genotype $index)
        (let $ntSet (start expr op var)
        (if (or (not (isMember $symbol $ntSet)) (== $index (len $genotype)))
            ($symbol $index)
            (let* 
                (
                    ($genotypeProb (selectByIndex $genotype $index))
                    (() (println! ("Next genotype prob" $genotypeProb)))
                    ($nextRule (getOption $grammarSpace $symbol $genotypeProb 0 0))
                    (() (println! ("Next rule" $nextRule)))
                    ; ($next (% (+ 1 $index) (len $genotype)))
                    ($next (+ 1 $index))
                )
                (ntByt $grammarSpace $nextRule $genotype $next)
            ))))

;; *********              ***********
;; ********* NEW APPROACH ***********
;; *********              ***********
;; *********              ***********

;; add a codon index in space 

! (add-atom &grammar (codon-index 0))

;; setGenotype  -- store the genotype in space for consistency of approach -- might be overkill
;;              $length -- length of the genotype tuple
;;              2 -- rounding

(= (setGenotype $grammarSpace $length)
    (let $genotype (createGenotype $length 2)
        (add-atom $grammarSpace (genotype $genotype))
    )
    
)

;; getCodon 
(= (getCodon $grammarSpace $dummy)
    (let*
        (
            ($index (match $grammarSpace (codon-index $x) $x))
            ($genotype (match $grammarSpace (genotype $g) $g))
            ($i (% $index 9)) ;; to enable wrapping around genotype 
            (() (println! ("Codon index " $index)))
        )
        (selectByIndex $genotype $i)
    )
)
(= (updateIndex $grammarSpace $rndm)
    (let*
        (
            ($index (match $grammarSpace (codon-index $x) $x))
            (() (println! ("index  " $index)))
            (() (remove-atom $grammarSpace (codon-index $index)))
            ($next (+ 1 $index))
        )
        (add-atom $grammarSpace (codon-index $next))
    )
)

;; generateIndividual expand a non-terminal into terminal only expressions

(= (generateIndividual $grammarSpace $start $random)
    (let*
        (
            ($codon (getCodon $grammarSpace ((py-atom random.random)))) ;; gets the corresponding codon in the geotype
            ($ntSet (start expr var op))
        )
        (if (and (isLiteral $start) (not (isMember $start $ntSet)))
            $start
            (if (== (get-metatype $start) Expression)
                (apply generateIndividual $grammarSpace $start 0 ((py-dot np random.uniform )0 1))
                (let*
                    (
                        ; (() (println! ("current rule" $start)))
                        ($nextRule (getOption $grammarSpace $start $codon 0 0))
                        (() (println! ("next codon" $codon)))
                        (() (println! ("next rule" $nextRule)))
                        (() (updateIndex $grammarSpace ((py-dot np random.uniform )0 1)))
                        (() (println! ("updated index " (match $grammarSpace (codon-index $x) $x))))
                    )
                    (generateIndividual $grammarSpace $nextRule ((py-atom random.random)))
                )))))
;; Helper function to apply the main function generateIndividual iteratively for derivation options like (expr op expr)
;; This will changed by using accumulator and append function -- faster execution

(= (apply $func $arg $expr $startIndex $dummy)
    (let* 
        (   ($len (length $expr))
            ; (() (println! ("length " $len)))
            ($e (selectByIndex $expr $startIndex))
            ($new ($func $arg $e ((py-atom random.random))))
            (() (println! ("new " $new)))
            ($updated (replaceByIndex $expr $startIndex $new))
            (() (println! ("updated " $updated)))
            ($current (+ 1 $startIndex))
        )
        (if (< $current $len)
            (apply $func $arg $updated $current ((py-dot np random.uniform )0 1))
            $updated
        )))
                ; ))
;; ************ RUN RUN ***************
;; ************ RUN RUN ***************
;; ************ RUN RUN ***************
;; ************ RUN RUN ***************

! (setGenotype &grammar 9)

! (generateIndividual &grammar start ((py-atom random.random)))

;; **************           **************
;; **************   TO      **************
;; **************     DO    **************
;; **************           **************

;; Optimization
;; generating population
;; ids for inidividuals and derivation option count for probablity update after scoring
;; scoring