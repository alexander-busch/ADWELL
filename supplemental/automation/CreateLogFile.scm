; Scheme function to create log file for Fluent usage
; http://www.eureka.im/835.html

	
; Often at the end of a FLUENT run we find it desirable if information for the run is available in the form of a log file which can be accessed later to give information about the case. Details such as version, date and time, cpu time, cortex time, number of nodes used and their details, models used etc are if written in a single file then it can serve as a quick and effective reference for the case. Attached scheme file does just that.

; This scheme file, if read at the end of the run provides a file named log-file.txt which contains the following information.
; 1. FLUENT version
; 2. Date and Time
; 3. Client Name
; 4. Cortex Name
; 5. Number of Compute Nodes
; 6. Names and IDs of the compute nodes used
; 7. Grid Information
; 8. Parallel Timer Information
; 9. Models used in FLUENT

; Here are the limitations of this function.
; 1. It works for FLUENT parallel version
; 2. It works only when UNIX based operating system is used as many system commands are being accessed.

; This file must be read once the iterations are complete for a FLUENT parallel version.



; This scm function creates a log file 'log-file.txt' where in report for mesh size, compute-nodes,
; CPU time, models, fluent version, date and time, client name, cortex name are included.
; This function works with parallel version of fluent.
; It must be run after iterations are complet to give the desired time log in the file
; This function works ONLY WITH UNIX/LINUX OPERATING SYSTEMS

(system "rm rabc-tf.txt")
(system "rm rabc-nf.trn rabc-n.trn rabc-n1.trn")
(system "rm rabc-mf.trn rabc-m1.trn rabc-m.trn")
(system "rm rabc-tf.trn rabc-t1.trn rabc-t.trn")
(system "rm rabc-mof rabc-mo1 rabc-mo")
(system "rm log-file.txt")
(system "cat rabc-tf.txt rabc-nf.trn rabc-mf.trn rabc-tf.trn rabc-mof > log-file.txt")
(append-file "rabc-tf.txt" `(lambda (port) (format port "~an","***************" )))
(append-file "rabc-tf.txt" `(lambda (port) (format port "~a","Fluent" )))
(append-file "rabc-tf.txt" `(lambda (port) (format port " ~an" (inquire-release))))
(append-file "rabc-tf.txt" `(lambda (port) (format port "~annn","***************" )))
(define da_ti "Date and Time")
(append-file "rabc-tf.txt" `(lambda (port) (format port "~a ", "Date and Time: " )))
(system "date >> rabc-tf.txt")
(append-file "rabc-tf.txt" `(lambda (port) (format port "n" )))
(append-file "rabc-tf.txt" `(lambda (port) (format port "~a ", "Client Name: " )))
(append-file "rabc-tf.txt" `(lambda (port) (format port " ~an" ,(cx-client-host))))
(append-file "rabc-tf.txt" `(lambda (port) (format port "~a", "Cortex Name: " )))
(append-file "rabc-tf.txt" `(lambda (port) (format port " ~annn" ,(cx-cortex-host))))
(append-file "rabc-tf.txt" `(lambda (port) (format port "~a ", "Compute Nodes = " )))
(append-file "rabc-tf.txt" `(lambda (port) (format port " ~annn" (compute-node-count))))
(system "cat rabc-tf.txt > log-file.txt")
(append-file "log-file.txt" `(lambda (port) (format port "~an ", "Compute Noes Used " )))
(define cmd-node "/file/start-transcript ")
(ti-menu-load-string (string-append cmd-node "rabc-n.trn" " y"))
(define cmd2-node "/parallel/show-connectivity")
(ti-menu-load-string (string-append cmd2-node " 1"))
(ti-menu-load-string "file stop-transcript")
(system "tail +3 rabc-n.trn >rabc-n1.trn")
(system "awk -F: '! /stop-transcript/' rabc-n1.trn > rabc-nf.trn")
(system "cat rabc-nf.trn >> log-file.txt")
(define cmd-mesh "/file/start-transcript ")
(ti-menu-load-string (string-append cmd-mesh "rabc-m.trn" " y"))
(define cmd2-mesh "/grid/size-info")
(ti-menu-load-string cmd2-mesh)
(ti-menu-load-string "file stop-transcript")
(system "tail -7 rabc-m.trn> rabc-m1.trn")
(system "awk -F: '! /stop-transcript/' rabc-m1.trn > rabc-mf.trn")
(system "cat rabc-mf.trn >> log-file.txt")
(append-file "log-file.txt" `(lambda (port) (format port "nn " )))
(define cmd-timer "/file/start-transcript ")
(ti-menu-load-string (string-append cmd-timer "rabc-t.trn" " y"))
(define cmd2-timer "/parallel/timer/print")
(append-file "rabc-t.trn" `(lambda (port) (format port "~ann","***************" )))
(ti-menu-load-string cmd2-timer)
(ti-menu-load-string "file stop-transcript")
(system "tail -18 rabc-t.trn > rabc-t1.trn")
(system "awk -F: '! /stop-transcript/' rabc-t1.trn > rabc-tf.trn")
(system "cat rabc-tf.trn >> log-file.txt")
(append-file "log-file.txt" `(lambda (port) (format port "nn " )))
(define cmd2-model-trash "/report/summary")
(ti-menu-load-string (string-append cmd2-model-trash " yes" " rabc-mo" " yes"))
(system "awk 'NR == 7, NR == 22' rabc-mo >rabc-mof")
(system "rm rabc-mo")
(system "cat rabc-mof >> log-file.txt")
(system "rm rabc-tf.txt")
(system "rm rabc-nf.trn rabc-n.trn rabc-n1.trn")
(system "rm rabc-mf.trn rabc-m1.trn rabc-m.trn")
(system "rm rabc-tf.trn rabc-t1.trn rabc-t.trn")
(system "rm rabc-mof rabc-mo1 rabc-mo")
