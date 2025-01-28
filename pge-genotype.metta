! (register-module! ./utilities)
! (import! &self utilities:general-helpers)

;; bind python's round function to rd
! (bind! np (py-atom numpy))
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
                ($first (rd ((py-dot np random.uniform)0 1) $prec))
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

! (bind! &grammar (new-space)) ;; for grammar use
! (bind! &geno (new-space)) ;; for genotype related

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
; ! (initP &grammar) ;; better call this after the length parameter of the setGenotype function has been set 

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
            
(= (getOption $grammarSpace $nt $genotypeProb $initialCumProb $index $dummy)
    (match $grammarSpace ($nt (($x $prob) $index)) 
            (if (<= $genotypeProb (+ $prob $initialCumProb)) 
                $x 
                (let* 
                    (
                        ($cumProb (+ $prob $initialCumProb))
                        ($next (+ $index 1))
                    )
                    (getOption $grammarSpace $nt $genotypeProb $cumProb $next $dummy)
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

;; add a codon index in space, &geno -- space for info related to genotype
;; genotype, codon index
! (add-atom &geno (codon-index 0))

;; setGenotype  -- store the genotype in space for consistency of approach -- might be overkill
;;              $length -- length of the genotype tuple
;;              2 -- rounding

(= (setGenotype $genoSpace $grammarSpace)
    (let* 
        (
            ($len (setLen $grammarSpace))
            ($genotype (createGenotype $len 2))
        )
        (add-atom $genoSpace (genotype $genotype))
    )
)
;; getRuleLength to use in the genotype creation -- should at least be long enough to allow for unique assignment for each rule in one pass

(= (setLen $grammarSpace)
    (sum (collapse (match $grammarSpace ($nt (($x) $index)) 1)))
)

; ! (setGenotype &geno &grammar)
; ! (match &geno (genotype $x) $x)
;; getCodon -- gets a codon in the genotype from 
(= (getCodon $genoSpace $dummy)
    (let*
        (
            ($index (match $genoSpace (codon-index $x) $x))
            ($genotype (match $genoSpace (genotype $g) $g))
            ($i (% $index 8)) ;; to enable wrapping around genotype 
            (() (println! ("Codon index " $index)))
        )
        (selectByIndex $genotype $i)
    )
)
;; updateIndex -- increments the codon-index value by 1 when called

(= (updateIndex $genoSpace $rndm)
    (let*
        (
            ($index (match $genoSpace (codon-index $x) $x))
            (() (println! ("index  " $index)))
            (() (remove-atom $genoSpace (codon-index $index)))
            ($next (+ 1 $index))
        )
        (add-atom $genoSpace (codon-index $next))
    )
)

;; generateIndividual expand a non-terminal into terminal only expressions

(= (generateIndividual $grammarSpace $genoSpace $start $random)
    (let*
        (
            ($codon (getCodon $genoSpace ((py-dot np random.uniform)0 1))) ;; gets the corresponding codon in the geotype
            ($ntSet (start expr var op))
        )
        (if (and (isLiteral $start) (not (isMember $start $ntSet)))
            $start
            (if (== (get-metatype $start) Expression)
                ; (iterate $func $arg1 $arg2 $expr $acc $dummy)
                (iterate generateIndividual $grammarSpace $genoSpace $start ((py-dot np random.uniform)0 1) ())
                ; (apply generateIndividual $grammarSpace $genoSpace $start 0 ((py-dot np random.uniform )0 1))
                (let*
                    (
                        ; (() (println! ("current rule" $start)))
                        ($nextRule (getOption $grammarSpace $start $codon 0 0 ((py-dot np random.uniform)0 1)))
                        (() (println! ("next codon " $codon)))
                        (() (println! ("next rule " $nextRule)))
                        (() (updateIndex $genoSpace ((py-dot np random.uniform)0 1)))
                        (() (println! ("updated index " (match $genoSpace (codon-index $x) $x))))
                    )
                    (generateIndividual $grammarSpace $genoSpace $nextRule ((py-dot np random.uniform)0 1))
                )))))
;; Helper function to apply the main function generateIndividual iteratively for derivation options like (expr op expr)
;; This will changed by using accumulator and append function -- faster execution

(= (apply $func $arg1 $arg2 $expr $startIndex $dummy)
    (let* 
        (   ($len (length $expr))
            ; (() (println! ("length " $len)))
            ($e (selectByIndex $expr $startIndex))
            ($new ($func $arg1 $arg2 $e ((py-dot np random.uniform))))
            (() (println! ("new " $new)))
            ($updated (replaceByIndex $expr $startIndex $new))
            (() (println! ("updated " $updated)))
            ($current (+ 1 $startIndex))
        )
        (if (< $current $len)
            (apply $func $arg1 $arg2 $updated $current ((py-dot np random.uniform)0 1))
            $updated
        )))

(= (iterate $func $arg1 $arg2 $expr $dummy $acc)
    (if (== $expr ())
        $acc
        (let* 
            (
                (($f $t) (decons $expr))
                ($r ($func $arg1 $arg2 $f ((py-dot np random.uniform)0 1)))
                ($c (append $acc $r))
                (() (println! ("current expression  " $c)))
            )
            (if (== $t ())
                $c
                (iterate $func $arg1 $arg2 $t ((py-dot np random.uniform)0 1) $c)
            ))))
                ; ))
;; ************ RUN RUN ***************
;; ************ RUN RUN ***************
;; ************ RUN RUN ***************
;; ************ RUN RUN ***************

! (setGenotype &geno &grammar)
! (initP &grammar)
! (generateIndividual &grammar &geno start ((py-dot np random.uniform)0 1))

;; **************           **************
;; **************   TO      **************
;; **************     DO    **************
;; **************           **************

;; Optimization
;; generating population
;; ids for inidividuals and derivation option count for probablity update after scoring
;; scoring

;; :::::::::::::::::::: Generate :::::::::::::::::
;; ::::::::::::::::::::::::::::: Population :::::::::::::::::

;; TO DO
;; generatePop function to call the generateIndividual function to generate popn of individuals
;; initialize the genotype of every individual
;; generate inidivdual -- set termination conditions for the individual generation (in generateIndividual)
;; store every fit individual in a space
;;      add id, genotype and derivation option count information together with the individual
;;      possible dict format -- ⤵️ 
;;                      (id 
;;                             (inidvidual $inidvidual)
;;                             (genotype $the_genotype)
;;                             (count 
;;                                   (
;;                                      (nt_0 ((rule_0 count_0))) (nt_1 ((rule_0 count_0) 
;;                                                                      (rule_1 count_1)
;;                                                                      (rule_2 count_2)
;;                                                                              :
;;                                                                                     ))
;;                                                            :
;;                                                            :
;;                                    )
;;                              )
;;                       )
;; Proposed Procedure: 
;;              main procedure genPop is called with spaces (grammar, popn space) and an accumulator with an id (=0) parameter assigned as (0)
;;              generate indivdual expression (counting of used derivation rules is done here)
;;              append the individual to (id) tuple -- (id (individual $individual))
;;              append the genotype used in the generation of the 
;;              individual $acc now --> (id (individual $individual) (genotype $genotype))
;;              get the updated count of all the derivation rules and append to acc
;;              acc now --> (id (individual $individual) (genotype $genotype) (count ()))
;;              add this tuple to space
;;              increment the id by 1, reset the derivation rule count and start the over untill popn is full 
;; Parameters  
;;          POP_SIZE
;;          Space to store individual expressions
;;          
;; derivation option counter for each valid expression
;;
;; setInitialCount -- a function to get and initialize count of derivation options to 0
;; should be able to reset counts from previous iteration
;; Note: the matching pattern may change depending on when the match is done -- 
;;          if this is done after the initial probability assignment the correct patten 
;;          is not ($nt (($rule) $x)) rather ($nt (($rule $prob) $x)) and vice versa
 
! (bind! &count (new-space))

;; set initial rule count -- 0 for each rule
(= (setInitialCount $grammarSpace $countSpace)
    (match $grammarSpace ($nt (($rule) $x)) 
        (add-atom $countSpace ($nt (($rule 0) $x)))
    ))

;; update rule count -- increment rule count by one when called

(= (updateRuleCount $countSpace $nt $rule)
    (let* 
        (
            (($count $index) (match $countSpace ($nt (($rule $c) $x)) ($c $x)))
            (() (remove-atom $countSpace ($nt (($rule $count) $index))))
            ($next (+ 1 $count))
        )
        (add-atom $countSpace ($nt (($rule $next) $index)))
    )
)
;; test setInitialCount

; ! (setInitialCount &grammar &count)
;; update a rule associated with an nt
; ! (updateRuleCount &count expr var)

;; test update
; ! (match &count (expr ((var $c)) $index) (var $c)) -- should be (var 1)

; ! (get-atoms &count)


(= (getRuleCount $countSpace $ntSet $acc)
    (let*
        (
            ($nt $t) (decons $ntSet)
            ($rulesCount (collapse (match $countSpace ($nt (($rule $count)) $index))))
            ($carry (append $acc ($nt $rulesCount)))
        )
        (if (== $t ())
            $carry
            (getRuleCount $countSpace $t $carry)
        )
    )
)

; ! (get-atoms &count) ; --> 4
;; genPop -- generates poulation and stores each individual to the popSpace
;; $id for the 
! (bind! &popn (new-space))

; (= (generateIndividual $grammarSpace $start $random)
; (= (getNt $grammarSpace)
; ! (setGenotype &grammar 9)
; ! (initP &grammar)

(= (genPop $popSize $grammarSpace $popSpace $id)

)