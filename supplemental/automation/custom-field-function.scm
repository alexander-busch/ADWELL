(custom-field-function/define
 '(((name solid-volumetric-flow-rate) (display "solid-vof * solid-x-velocity") (syntax-tree ("*" "solid-vof" "solid-x-velocity")) (code (field-* (field-load "solid-vof") (field-load "solid-x-velocity"))))
   ))