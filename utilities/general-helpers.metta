;; Lenght of a tuple

(= (len $tuple)
    (if (== $tuple ())
        0
        (let $tail (cdr-atom $tuple)
            (+ 1 (len $tail))
        )))
;; non-deterministic definitio of length

(= (length $expr)
    (sum (collapse (let $a (superpose $expr) 1)))
)
(: isSymbol (-> Expression Bool))
(= (isSymbol $expr) (== (get-metatype $expr) Symbol))

(: isUnit (-> Expression Bool))
(= (isUnit $expr) (== $expr ()))

(: isSingleAtom (-> Expression Bool))
(= (isSingleAtom $expr) (== (cdr-atom $expr) ()))

;; Function to replace a given atom from a space.
(: update-atom (-> Grounded Atom Atom (->)))
(= (update-atom $space $oatom $natom) (let $_ (remove-atom $space $oatom) (add-atom $space $natom)))

;; A function to return head and tail of an atom in a tuple.
(= (decons $atom)
   (if (not (or (isSymbol $atom) (isUnit $atom)))
       ((car-atom $atom) (cdr-atom $atom))
       (Error Unit Or Symbol cannot be deconstructed)))

;; Fold a tuple from right to left
(: foldr (-> (-> $a $b $b) $b $c $d))
(= (foldr $f $i $xs)
   (if (== $xs ())
       $i
       (let* (($h (car-atom $xs))
              ($t (cdr-atom $xs))
              ($ft (foldr $f $i $t)))
         ($f $h $ft))))

;; Add numbers in an atom list
(= (sum $list) (foldr + 0 $list))

;; Select atom by index

(= (selectByIndex $expr $index)
    (if (== $expr ())
        (Error (Index out of range) ())
        (if (== $index 0)
            (car-atom $expr)
            (let $tail (cdr-atom $expr)
                (selectByIndex $tail (- $index 1))
            ))))

;; deletes atom by index

(= (deleteAtom $expr $index)
    (if (> $index (- (len $expr) 1))
        (Error (Index out of range) ())
        (if (== $index 0)
            (cdr-atom $expr)
            (let* 
                (
                    (($first $tail) (decons $expr))
                    ($carry (deleteAtom $tail (- $index 1)))
                )
                (cons-atom $first $carry)
            ))))

;; isLiteral -- determines if the given atom is a literal or not
(= (isLiteral $a)
    (if (or (== (get-metatype $a) Grounded) (== (get-metatype $a) Symbol))
        True
        (if (== (car-atom $a) NOT)
            True
            False )))

; a function which checks if an element is member of a tuple.
(: isMember (-> $a $aa Bool))
(= (isMember $x $list)
    (not (== (collapse (intersection (superpose ($x)) (superpose $list))) ()))
)

;; A function that behaves like a do while loop.
;; It executes the functions in the given list and returns the last result.
(: until (-> (-> $a Bool) (-> $a $a) $a $a))
(= (until $predicate $function $x) 
    (if ($predicate $x) 
        $x 
        ($y ($function $x) (until $predicate $function $y))
    ))

;; a function which returns index of first occurence of an atom in a tuple which satisfies a condition
;; 0 indexed counting

(= (findX $predicate $expr)
    (if (== $expr ())
        (Error (No such element!) ())
        (let ($f $t) (decons $expr)
            (if ($predicate $f)
                0
                (let $r (findX $predicate $t)
                    (if (== (get-type $r) Number)
                        (+ 1 $r)
                        $r
                    ))))))

; ! (findX isLiteral ((expr) expr (expr) (expr) expr (expr)))

;; replaceByIndex replaces an atom at a specific index position with $val -- 0 indexed counting is implied

(: replaceByIndex (-> Expression Number $t Expression))
(= (replaceByIndex $expr $index $new)
    (if (== $expr ())
        $expr
        (let*
            (
                ($first (car-atom $expr))
                ($tail (cdr-atom $expr))
            )
            (if (== $index 0)
                (if (== $tail ())
                    ($new)
                    (cons-atom $new $tail)
                )
                (let $c (replaceByIndex $tail (- $index 1) $new)
                    (cons-atom $first $c)
                )))))

;; apply a function by selecting and replacing using index
;; map doesn't allow for keeping the index of genotype 

(= (iterativeApply $func $y $expr $index)
    (let $len (len $expr)
        (if (< $index $len)
            (let* 
                (
                    ($e (selectByIndex $expr $index))
                    ($r ($func $e $y))
                    ($update (replaceByIndex $expr $index $r))
                    ($next (+ 1 $index))
                )
                (iterativeApply $func $y $update $next)
            )
            $expr
        )
    )
)

(= (map $func $dummy $expr)
    (if (isLiteral $expr)
        ($func $dummy $expr)
        (let*
            (
                (($f $t) (decons $expr))
                ($c (map $func $dummy $f))
            )
            (if (== $t ())
                (cons-atom $c $t)
                (let $d (map $func $dummy $t)
                    (cons-atom $c $d)
                )))))

;; append atom to the end of an expression

(= (append $acc $x)(collapse (union (superpose $acc) $x)))

;; min finds minimum of two numbers
(: min (-> Number Number Number))
(= (min $x $y) (if (<= $x $y) $x $y))

;; slice -- gets a range of values within the specified start and end indices values
;;          $expr -- tuple to be sliced
;;          $s -- start index of the slicer
;;          $e -- end index of the slicer
;;          $acc -- an empty accumulator -- ** an implementation without an accumulator should be possible **
;;          start and end indices check could be added to make sure the indices are not out of bounds

(= (slice $expr $s $e $acc)
    (let ($f $t) (decons $expr)
        (if (== $s 0)
            (let $c (append $acc $f)
                (if (== $e 0)
                    $c
                    (slice $t 0 (- $e 1) $c)
                )
            )
            (slice $t (- $s 1) (- $e 1) $acc)
        )))

(= (max $x $y) (if (>= $x $y) $x $y))

;; maxIter -- finds maximum from a tuple of numerals

(= (maxIter $expr)
    (let ($f $t) (decons $expr)
        (if (== $t ())
            $f
            (max $f (maxIter $t))
        )
    )
)


; Concatenates two tuples $xs and $ys e.g (concatTuple (1 2 3) (4 5)) -> (1 2 3 4 5)

(= (concatTuple $xs $ys ) (collapse (union (superpose $xs) (superpose $ys))))

;; This function slices a tuples at the given index and returns a tuple of the sliced part and the remaining tuple.
;;              Example: (slice (1 2 3 4 5) 2 ()) -> ((1 2 3) (4 5))
(: slice (-> Expression Number Expression Expression))
(= (slice $lst $index $acc)
    (if (or (== $lst ()) (< $index 0))
        $acc
        (let*
            (
                ($head (car-atom $lst))
                ($tail (cdr-atom $lst))
                ($newAcc (append $acc $head))
                
            )
            (if (== $index 0)
                ($newAcc $tail)
                (slice $tail (- $index 1) $newAcc)
            )
                
        )
    )
)
