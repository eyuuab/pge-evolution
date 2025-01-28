! (register-module! ./utilities)
! (import! &self utilities:general-helpers)

;; bind python's round function to rd
! (bind! rd (py-atom round))

;; space for grammar rules
! (bind! &grammar (new-space))

;; add the grammar to &grammar space
! (add-atom &grammar (start ((expr) 0)))
! (add-atom &grammar (expr (((expr op expr)) 0)))
! (add-atom &grammar (expr ((var) 1)))
! (add-atom &grammar (op (("+") 0)))
! (add-atom &grammar (op (("*") 1)))
! (add-atom &grammar (op (("-") 2)))
! (add-atom &grammar (var ((x) 0)))
! (add-atom &grammar (var ((1.0) 1)))

;; &count -- space for derivation rules usage count 
! (bind! &count (new-space))

;; setInitialCount -- initialize rule count to 0
(= (setInitialCount $grammarSpace $countSpace)
    (match $grammarSpace ($nt (($rule) $index)) 
        (add-atom $countSpace ($nt (($rule 0) $index)))
    ))

;; updateRuleCount -- increments derivation rule count of a given non-terminal by one when called
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
! (setInitialCount &grammar &count)

;; updateRuleCount -- sample call
! (updateRuleCount &count expr var)

;; getNt -- returns tuple of unique non-terminals of the grammar 
(= (getNt $grammarSpace)
    (collapse (unique (match &grammar ($nt $x) $nt)))
)

;; getRuleCount -- returns count of derivation rules as tuple of tuples
;; Example -- type of the return tuple according to the given grammar
;;            ((start (expr 0)) (expr ((expr op expr) 1) (var 0)) ...)
(= (getRuleCount $countSpace $ntSet $acc)
    (let*
        (
            (($nt $t) (decons $ntSet))
            ($rulesCount (collapse (match $countSpace ($nt (($rule $count) $index)) ($rule $count))))
            ($carry (append $acc ($nt $rulesCount)))
        )
        (if (== $t ())
            $carry
            (getRuleCount $countSpace $t $carry)
        )
    )
)

;; updateProbablity -- updates probability of a derivation rules based on count information from &count space 
(= (updateProbablity $grammarSpace $countSpace $lambda)
    (let* 
        (
            ($nts (getNt $grammarSpace))
            ($nt (superpose $nts))
            ($c (getRuleCount $countSpace $nts ()))
            ($sums (getCountSum $c $nt))
        )
        (match $grammarSpace ($nt (($rule $prob) $index)) 
            (match $countSpace ($nt (($rule $count) $index))
                (if (== $count 0.0)
                    (let*
                        (
                            ($dec (* $lambda $prob))
                            ($newProb (rd (- $prob $dec) 3))
                            (() (remove-atom $grammarSpace ($nt (($rule $prob) $index))))
                        )
                        (add-atom $grammarSpace ($nt (($rule $newProb) $index)))
                    )
                    (let*
                        (   
                            ($temp (* $count $lambda))
                            ($inc (/ $temp $sums))
                            ($incP (+ $prob $inc))
                            ($newProb (rd (min $incP 1) 3))
                            (() (remove-atom $grammarSpace ($nt (($rule $prob) $index))))
                        )
                        (add-atom $grammarSpace ($nt (($rule $newProb) $index)))
                    )
                )))))

;; getCountSum -- returns sum of counts of derivation ruels for each non-terminal
(= (getCountSum $allCount $nt)
    (let $cs (superpose $allCount)
            (unify $cs ($nt $x) 
                    (sum (collapse (let $rule (superpose $x)
                        (unify $rule ($y $count)
                            $count
                            (empty)
                        ))))
                    (empty)
            )))
;; :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;; ::::::::::::::  ASSIGN                      :::::::::::::
;; ::::::::::::::        INITIAL               :::::::::::::
;; ::::::::::::::               PROBABILITIES  :::::::::::::
;; :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
;;
(= (updateRule $grammarSpace $nt $len)
    (let () (match $grammarSpace ($nt (($x) $n)) (add-atom $grammarSpace ($nt (($x (rd (/ 1 $len) 2)) $n))))
    
        (match $grammarSpace ($nt (($x) $n)) (remove-atom $grammarSpace ($nt (($x) $n))))
    )
)
;; ************ initialize initial probabality for grammar derivation options ***********

(= (initP $grammarSpace)
    (let $ntSet (getNt $grammarSpace)
        (assignInitialProbs $grammarSpace $ntSet)
    )
)
! (initP &grammar)

;; ::::::::::::::::: END OF INITIAL PROBABILITY ASSIGNMENT :::::::::::::::::

;; :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;; ::::::::::::::  RUN                         :::::::::::::
;; ::::::::::::::        UPDATE                :::::::::::::
;; ::::::::::::::               PROBABILITIES  :::::::::::::
;; :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! (updateProbablity &grammar &count 0.05)

;; CHECK THE UPDATE 
;; unused rules will have their probablity values decreased while used rules' probablity values will be increased
! (get-atoms &grammar)

;; normalize -- adjusts probabilities of derivation rules sum to be exactly one
;;        chucks if sum of probablities of all rules is greater than one
;;        if so, subtract from all derivation rules the same amount given by the following formula  

(= (normalizeP $grammarSpace)
    (let*
        (
            ($ntSet (getNt $grammarSpace))
            ($nt (superpose $ntSet))
            ($probSums (sum (collapse (match $grammarSpace ($nt (($rule $prob) $index)) $prob))))
            ($j (sum (collapse (match $grammarSpace ($nt (($rule $prob) $index)) 1))))
            ($diff (rd (/ (- 1 $probSums) $j) 3))
            ((($oldProb $newProb) ($rule $index)) (match $grammarSpace ($nt (($rule $prob) $index)) (($prob (rd (+ $prob $diff) 3)) ($rule $index))))
            ; ($oldProb (match $grammarSpace ($nt (($rule $prob) $index))))
            (() (println! ("new P  " $newProb)))
            (() (remove-atom $grammarSpace ($nt (($rule $oldProb) $index))))
        )
        
        (add-atom $grammarSpace ($nt (($rule $newProb) $index)))
    ))
;; ************** RUN **************
;; ************** NORMALIZE **************
;; ************** PROBABILITY **************

! (normalizeP &grammar)

;; ************** Check normalized probabilities ************
;; sum of probabilities of rules of a non-terminal must equal 1
! (get-atoms &grammar)



