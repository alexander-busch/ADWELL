;;=====================================================;;; 
;;; This Scheme function is computing 
;;; pressure gradient of a translational periodic 
;;; 
;;; Source http://www.eureka.im/295.html
;;; 
;;; Usage 

;;; 1) Read the case and data file, 
;;; 
;;; 2) File>Read>Scheme .. this file (monitor-periodic.scm) 
;;; In order to load this file automatically 
;;; Create a .fluent file in your home directory in which you type the following line: 
;;; (load "~/mypath/monitor-periodic.scm") 

;;; 3) Set a monitor command in GUI 
;;; Solve->Execute Commands 
;;; with the following content 
;;; 
;;; (write-my-deltap filename) 
;;; 
;;; filename: a string containing the name of your output file 
;;; 
;;; Example: 
;;; (write-my-deltap "monitor1.out" ) 
;;; 

;;; 4) Set the frequency of executions - Every "n" Time steps or iterations 

;;; 5) Run 
;;; 
;;; Do not change below this line 
;;; ========================================================== 
(define (write-my-deltap filename) 

(let ((uport)(niter (%iterate 0))) 

(if (file-exists? filename) 

(set! uport (open-file filename "a+")) 
(let* ((fn filename)) 
(set! uport (open-file fn "w+")) 
(format uport "Pressure gradient as a function of iteration") 
(newline uport) 
) 
) 
(format uport " ~a" niter) 
(format uport " t ~a" (rpgetvar 'periodic/pressure-derivative)) 
(newline uport) 
(flush-output-port uport) 
(close-output-port uport) 
) 
) 