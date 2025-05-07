
set pi 3.1415926

proc evaluate_usage { } {

  vmdcon -info "Usage: evaluate <psffile> <dcdfile> -ref <reference> \[options...\]" 
  vmdcon -info "Calculate hydration layer structure"
  vmdcon -info "Required Parameters:"
  vmdcon -info "  <psffile>       Structure file"
  vmdcon -info "  <dcdfile>       Trajectory file"
  vmdcon -info "  -ref <sel>      Reference atom selection (required)"
  vmdcon -info "Options:"
  vmdcon -info "  -first <int>    First frame (default: 0)"
  vmdcon -info "  -last <int>     Last frame (default: -1)"
  vmdcon -info "  -n <int>        Number of segments (default: 4)"
  vmdcon -info "  -s <int>        Stride value (default: 1)"
  vmdcon -info "  -cal <sel>      Calculation atoms (default: type OT)"
  vmdcon -info "  -r <float>      Hydration radius (default: 3.5)"
  vmdcon -info "  -RDF <0/1>      Enable RDF calculation (default: 0)"
  vmdcon -info "  -ODF <0/1>      Enable ODF calculation (default: 0)"
  vmdcon -info "  -HB  <0/1>      Enable HB  analysis    (default: 0)"
  vmdcon -info "  -MSD <0/1>      Enable MSD calculation (default: 0)"
  vmdcon -info "  -MRT <0/1>      Enable MRT calculation (default: 0)"
  error ""

}

proc evaluate {args} {

  global errorInfo errorCode

  set oldcontext [psfcontext new]
  set errflag [catch { eval evaluate_core $args } errMsg]
  set savedInfo $errorInfo
  set savedCode $errorCode
  psfcontext $oldcontext delete
  if {$errflag} { error $errMsg $savedInfo $savedCode }
	
}

proc evaluate_core {args} {

  set vmdexec vmd

  # Set some defaults

    # Print usage information if no arguments are given

    if {[llength $args] < 2} { evaluate_usage }

    # Confirm the existence of psffile and dcdfile

    lassign $args psffile dcdfile
    set args [lrange $args 2 end]

    if {![file exists $psffile]} { error "psffile not found: $psffile" }
    if {![file exists $dcdfile]} { error "dcdfile not found: $dcdfile" }

    # Check for even number of args

    if {[llength $args] % 2 != 0} { evaluate_usage }

    # Specify initialization options

    array set options {

      -first 0    -last -1     -n 4      -s 1
      -cal "type OT"           -r 3.5
      -RDF 0      -ODF 0       -HB 0     -MSD 0     -MRT 0

    }

    array set user_params $args

    set required_params {-ref}

    foreach param $required_params {

      if {![info exists user_params($param)]} { error "Missing required parameter: $param" }

    }

    array set options $args

    proc ::validate_int {val} {string is integer -strict $val}
    proc ::validate_float {val} {string is double -strict $val}
    proc ::validate_switch {val} {expr {$val in {0 1}}}

    foreach {key validator} {

      -first  validate_int    -last   validate_int
      -n      validate_int    -s      validate_int
      -r      validate_float
      -RDF    validate_switch -ODF    validate_switch
      -HB     validate_switch -MSD    validate_switch
      -MRT    validate_switch
    
    } {

      if {![{*}$validator $options($key)]} { error "Invalid value for $key: '$options($key)'" }

    }

  # Start counting

    # Split DCD File

    set step [splitDCD $psffile $dcdfile $options(-first) $options(-last) $options(-n) $options(-s)]

    # Dump scripts

    dumpscripts $psffile $step $options(-n) "$options(-ref)" "$options(-cal)" $options(-r) $options(-RDF) $options(-ODF) $options(-HB) $options(-MSD) $options(-MRT)

    # Parallel run

    parallelrun $options(-n) $vmdexec $step

  # Check and integrate

    # Make check and integration

    foreach analysis {RDF ODF HB MSD MRT} {

      if {$options(-$analysis)} {

        check $options(-n) "Data_$analysis"
        integrate_$analysis $options(-n)
      }

    }

  # Delete temporary files

    for {set i 0} {$i < $options(-n)} {incr i} {

      file delete $i.dcd
      file delete $i.tcl

    }

}


proc splitDCD {psffile dcdfile first last n stride} {

  if {$last != -1 && $last < $first} {error "Invalid frame range: last($last) < first($first)"}

  mol load psf $psffile
  mol addfile $dcdfile type dcd first $first last $last step $stride waitfor all
  set nframe [molinfo top get numframes]

  if {[molinfo top get numframes] == 0} {error "Couldn't load psf/pdb files!"}

  set id [molinfo top get id]

  if {$nframe % $n != 0} {error "The number of frames in each segment must be an integer"}

  set step [expr $nframe / $n]

  vmdcon -info "Splitting DCD into $n segments..."

  for { set i 0 } { $i < $n } { incr i 1 } {

    set nb [expr $i * $step]
    set ne [expr ($i + 1) * $step - 1]
    vmdcon -info "$i.dcd: from frame $nb to frame $ne"
    animate write dcd "$i.dcd" beg $nb end $ne $id

  }

  vmdcon -info "Note: above files will be deleted later."
  mol delete all
  return $step

}

proc dumpscripts {psffile step n refatm calatm r RDF ODF HB MSD MRT} {

  vmdcon -info "writing scripts for each trajectory ..."

  for {set i 0} {$i < $n} {incr i} {

    set initframe [expr $i * $step]
    set script [open $i.tcl w]

    puts $script "mol load psf $psffile dcd $i.dcd"
    puts $script "source Evaluator.tcl"

    if { $RDF } {puts $script "measure_RDF Data_RDF_$i \"$refatm\" \"$calatm\""}
    if { $ODF } {puts $script "measure_ODF Data_ODF_$i \"$refatm\" \"$calatm\" $r"}
    if { $HB } {puts $script "measure_HB Data_HB_$i \"$refatm\" \"$calatm\" $r"}
    if { $MSD } {puts $script "measure_MSD Data_MSD_$i \"$refatm\" \"$calatm\" $r"}
    if { $MRT } {puts $script "measure_MRT Data_MRT_$i \"$refatm\" \"$calatm\" $r"}

    puts $script "mol delete all"
    puts $script "exit"
    puts $script "   "

    close $script
    vmdcon -info "script $i.tcl has been written."

  }

}

proc parallelrun {n vmdexec step} {

  puts "Running tasks ..."
  
  # First (n - 1) tasks will be run in background.
  
  for { set i 0 } { $i < [expr $n - 1] } { incr i 1 } {

    exec vmd -dispdev none -e $i.tcl > /dev/null 2>> /dev/null &

  }

  # The last task will be run in foreground. When this task
  # terminaled, all tasks should have been terminaled.

  set lasttask [expr $n - 1]
  exec $vmdexec -dispdev none -e $lasttask.tcl > /dev/null 2>> /dev/null

  # Wait

  sleep [expr $n * 1.5]
  puts "#--------------------------------------"

}

proc check {n fileid} {

    set flag 0

    for {set i 0} {$i < $n} {incr i} {

      if {[file exists $fileid\_$i.dat]} {

        set flag [expr $flag + 1]

      }

    }

    if {$flag == $n} {

      puts "$fileid tasks finished."

    } else {

      error "$fileid tasks unfinished!!!"

    }

}

proc measure_RDF {outputname refatm {calatm "type OT"} {dr "0.1"} {rmax "10.0"}} {

  global pi

  set fileid [open $outputname\.dat w]

  set nframe [molinfo top get numframes]

  set refsel [atomselect top "$refatm"]
  set numref [$refsel num]
  set calsel [atomselect top "$calatm"]
  set numcal [$calsel num]

  set numbin [expr int(double($rmax)/$dr)]

  set l [molinfo top get a]
  set w [molinfo top get b]
  set h [molinfo top get c]
  set V [expr $l * $w * $h]
  set avgp_cal [expr double($numcal) / $V]

  array set num {}
  array unset num *

  for {set j 0} {$j < $numbin} {incr j} {

    set rmid [expr ($j + 0.5) * $dr]
    set num($rmid) 0

  }

  for {set i 0} {$i < $nframe} {incr i} {

    # puts -nonewline [format "      Now is calculating %i_th frames.....\n" $i]

    set rtmp 0

    for {set j 0} {$j < $numbin} {incr j} {

      set rmid [expr ($j + 0.5) * $dr]

      foreach k [$refsel get index] {

        set sel1 [atomselect top "$calatm && pbwithin $rtmp of index $k"]
        set sel2 [atomselect top "$calatm && pbwithin [expr $rtmp + $dr] of index $k"]

        $sel1 frame $i
        $sel2 frame $i

        $sel1 update
        $sel2 update

        set num($rmid) [expr [$sel2 num] - [$sel1 num] + $num($rmid)]

        $sel1 delete
        $sel2 delete

      }

      set rtmp [expr ($j + 1) * $dr]

    }

  }

  for {set j 0} {$j < $numbin} {incr j} {

    set rtmp [expr $j * $dr]
    set rmid [expr ($j + 0.5) * $dr]
    set dn   [expr double($num($rmid)) / (double($nframe) * double($numref))]
    set dV   [expr 4.0/3.0 * $pi * (3.0 * $rtmp**2 * $dr + 3.0 * $rtmp * ($dr)**2 + ($dr)**3)]
    set dp   [expr $dn / $dV]
    puts $fileid "$rmid [expr $dp / $avgp_cal]"

  }

  $refsel delete
  $calsel delete
  close $fileid
  
}

proc integrate_RDF {n} {

  set first_column [list]
  set all_data [list]

  set outfile [open "Data_RDF.dat" w]

  for {set i 0} {$i < $n} {incr i} {

    set fd [open Data_RDF_$i\.dat r]
    set lines [list]

    while {[gets $fd line] >= 0} {

      set columns [regexp -all -inline {\S+} $line]

      if {$i == 0} {

        lappend first_column [lindex $columns 0]

      }

      lappend lines [lindex $columns 1]

    }

    close $fd
    lappend all_data $lines

  }

  set num_rows [llength $first_column]

  for {set i 0} {$i < $num_rows} {incr i} {

    set sum 0.0
    for {set j 0} {$j < $n} {incr j} {

      set yval [lindex [lindex $all_data $j] $i]
      set sum [expr {$sum + $yval}]

    }

    set avg [expr {$sum / $n}]
    puts $outfile "[lindex $first_column $i] $avg"

  }

  close $outfile
}

proc measure_ODF {outputname refatm {calatm "type OT"} {rmax "3.5"}} {

  global pi

  set fileid [open $outputname\.dat w]

  set nframe [molinfo top get numframes]

  set refsel [atomselect top "$refatm"]

  array set frequency {}
  array unset frequency *

  set totalnum 0

  for {set i 1} {$i <= 180} {incr i} {

    set frequency($i) 0
    
  }

  for {set i 0} {$i <= $nframe} {incr i} {

    # puts -nonewline [format "      Now is calculating %i_th frames.....\n" $i]

    $refsel frame $i
    $refsel update

    foreach j [$refsel get index] {

      set sel_water [atomselect top "$calatm && pbwithin $rmax of index $j"]
      $sel_water frame $i
      $sel_water update

      foreach k [$sel_water get index] {

        set sel_OT [atomselect top "index $k"]
        set sel_HT [atomselect top "index [lindex [lindex [$sel_OT getbonds] 0] 0] || index [lindex [lindex [$sel_OT getbonds] 0] 1]"]
        set sel_tp [atomselect top "index $j"]

        $sel_OT frame $i
        $sel_HT frame $i
        $sel_tp frame $i

        $sel_OT update
        $sel_HT update
        $sel_tp update

        set VecWat [vecsub [measure center $sel_HT] [measure center $sel_OT]]
        set Vectmp [vecsub [measure center $sel_OT] [measure center $sel_tp]]

        set Cos    [expr [vecdot $VecWat $Vectmp]/([veclength $VecWat]*[veclength $Vectmp])]
        set Angtmp [expr acos($Cos)*180/$pi]
        set Ang    [expr int(ceil($Angtmp))]

        incr frequency($Ang)
        incr totalnum

        $sel_OT delete
        $sel_HT delete
        $sel_tp delete

      }

      $sel_water delete

    }

  }

  puts $fileid "$totalnum"

  for {set i 1} {$i <= 180} {incr i} {

    puts $fileid "$frequency($i)"
    
  }

  $refsel delete
  close $fileid

}

proc integrate_ODF {n} {

  array set data {}
  array unset data *

  set outfile [open "Data_ODF.dat" w]

  for {set i 0} {$i < $n} {incr i} {

    set fd [open Data_ODF_$i\.dat r]
    set row 0

    while {[gets $fd line] >= 0} {

      set columns [regexp -all -inline {\S+} $line]

      lappend data($row) [lindex $columns 0]
      incr row

    }

    close $fd
        
  }

  set totalnum 0

  for {set i 0} {$i < 181} {incr i} {

    set value 0

    for {set j 0} {$j < $n} {incr j} {

      if {$i == 0} {

        set totalnum [expr [lindex $data($i) $j] + $totalnum]

      } else {

        set value [expr [lindex $data($i) $j] + $value]

      }
          
    }

    if {$i == 0} {

      puts $outfile "$totalnum"

    } else {

      puts $outfile "[expr double($value) / double($totalnum) * 100]"

    }
        
  }

  close $outfile
  
}

proc measure_HB {outputname refatm {calatm "type OT"} {rmax "3.5"} {amax "30"}} {

  global pi

  set fileid [open $outputname\.dat w]

  set nframe [molinfo top get numframes]

  set refsel [atomselect top "$refatm"]

  puts $fileid "$nframe [$refsel num]"

  set radians [expr {$amax * $pi / 180.0}]
  set cos_value [expr {cos($radians)}]

  for {set i 0} {$i < $nframe} {incr i} {

    # puts -nonewline [format "      Now is calculating %i_th frames.....\n" $i]

    $refsel frame $i
    $refsel update

    # now is calculating the ref atoms as HB acceptor

    foreach j [$refsel get index] {

      set sel_water [atomselect top "$calatm && pbwithin $rmax of index $j"]
      $sel_water frame $i
      $sel_water update

      foreach k [$sel_water get index] {

        set sel_OT  [atomselect top "index $k"]
        set sel_HT1 [atomselect top "index [lindex [lindex [$sel_OT getbonds] 0] 0]"]
        set sel_HT2 [atomselect top "index [lindex [lindex [$sel_OT getbonds] 0] 1]"]
        set sel_tp  [atomselect top "index $j"]

        $sel_OT  frame $i
        $sel_HT1 frame $i
        $sel_HT2 frame $i
        $sel_tp  frame $i

        $sel_OT  update
        $sel_HT1 update
        $sel_HT2 update
        $sel_tp  update

        set Vec1 [vecsub [measure center $sel_HT1] [measure center $sel_OT]]
        set Vec2 [vecsub [measure center $sel_HT2] [measure center $sel_OT]]
        set Vec3 [vecsub [measure center $sel_tp]  [measure center $sel_OT]]

        set Cos1 [expr [vecdot $Vec1 $Vec3]/([veclength $Vec1]*[veclength $Vec3])]
        set Cos2 [expr [vecdot $Vec2 $Vec3]/([veclength $Vec2]*[veclength $Vec3])]

        if {$Cos1 > $cos_value} {puts $fileid "1 [veclength $Vec3] [expr acos($Cos1)*180/$pi]"}
        if {$Cos2 > $cos_value} {puts $fileid "1 [veclength $Vec3] [expr acos($Cos2)*180/$pi]"}

        $sel_OT  delete
        $sel_HT1 delete
        $sel_HT2 delete
        $sel_tp  delete

      }

      $sel_water delete

    }

    # now is calculating the ref atoms as HB donor

    foreach j [$refsel get index] {

      set sel_tp [atomselect top "index $j"]

      $sel_tp frame $i
      $sel_tp update

      set list  [lindex [$sel_tp getbonds] 0]
      set llist [llength $list]

      for {set m 0} {$m < $llist} {incr m} {

        set sel [atomselect top "index [lindex $list $m]"]

        $sel frame $i
        $sel update

        if {[expr floor([$sel get mass])] == 1} {

          set sel_water [atomselect top "$calatm && pbwithin $rmax of index $j"]
          $sel_water frame $i
          $sel_water update

          foreach k [$sel_water get index] {

            set sel_OT [atomselect top "index $k"]

            $sel_OT frame $i
            $sel_OT update

            set Vec1 [vecsub [measure center $sel] [measure center $sel_tp]]
            set Vec2 [vecsub [measure center $sel_OT] [measure center $sel_tp]]

            set Cos1 [expr [vecdot $Vec1 $Vec2]/([veclength $Vec1]*[veclength $Vec2])]

            if {$Cos1 > $cos_value} {puts $fileid "0 [veclength $Vec2] [expr acos($Cos1)*180/$pi]"}

            $sel_OT delete
            
          }

          $sel_water delete
         
        }

        $sel delete
        
      }

      $sel_tp delete
      
    }

  }

  $refsel delete
  close $fileid

}

proc integrate_HB {n {rmax "3.5"} {dr "0.1"} {amax "30"} {da "1"}} {

  set outfile [open "Data_HB.dat" w]

  set nframe 0
  set numref 0
  set count(0) 0
  set count(1) 0

  set num_dist_bins  [expr {int(ceil($rmax / $dr))}]
  set num_angle_bins [expr {int(ceil($amax / $da))}]

  array set dist_bins  {0 {} 1 {}} 
  array set angle_bins {0 {} 1 {}}

  for {set t 0} {$t <= 1} {incr t} {

    for {set i 0} {$i < $num_dist_bins} {incr i} {

      lappend dist_bins($t) 0

    }

    for {set i 0} {$i < $num_angle_bins} {incr i} {

      lappend angle_bins($t) 0

    }

  }
  
  for {set i 0} {$i < $n} {incr i} {

    set fd [open Data_HB_$i\.dat r]

    gets $fd line

    set header [regexp -all -inline {\S+} $line]
    set frames [lindex $header 0]
    set atoms [lindex $header 1]

    if {$i == 0} {set numref $atoms}

    incr nframe $frames

    while {[gets $fd line] >= 0} {

      set values [regexp -all -inline {\S+} $line]
      set type [lindex $values 0]
      set distance [lindex $values 1]
      set angle [lindex $values 2]

      incr count($type)

      set dist_bin [expr {int(floor($distance / $dr))}]

      set angle_bin [expr {int(floor($angle / $da))}]

      lset dist_bins($type) $dist_bin [expr {[lindex $dist_bins($type) $dist_bin] + 1}]
      lset angle_bins($type) $angle_bin [expr {[lindex $angle_bins($type) $angle_bin] + 1}]

    }

    close $fd

  }

  set norm_factor [expr $nframe * $numref]
  set hb_acceptor [expr {$norm_factor ? double($count(1)) / $norm_factor : 0.0}]
  set hb_donor    [expr {$norm_factor ? double($count(0)) / $norm_factor : 0.0}]

  puts $outfile [format "%.4f %.4f" $hb_acceptor $hb_donor]

  puts $outfile "# Distance Distribution"

  for {set i 0} {$i < $num_dist_bins} {incr i} {

    set prob_acceptor [expr {($count(1) ? [lindex $dist_bins(1) $i] / double($count(1)) : 0.0) * 100}]
    set prob_donor    [expr {($count(0) ? [lindex $dist_bins(0) $i] / double($count(0)) : 0.0) * 100}]

    puts $outfile [format "%.3f %.3f" $prob_acceptor $prob_donor]

  }

  puts $outfile "# Angle Distribution"

  for {set i 0} {$i < $num_angle_bins} {incr i} {

    set prob_acceptor [expr {($count(1) ? [lindex $angle_bins(1) $i] / double($count(1)) : 0.0) * 100}]
    set prob_donor    [expr {($count(0) ? [lindex $angle_bins(0) $i] / double($count(0)) : 0.0) * 100}]

    puts $outfile [format "%.3f %.3f" $prob_acceptor $prob_donor]

  }

  close $outfile
  
}


proc measure_MSD {outputname refatm {calatm "type OT"} {rmax "3.5"}} {

  set fileid [open $outputname\.dat w]

  set nframe [molinfo top get numframes]

  set refsel [atomselect top "$refatm"]

  for {set i 0} {$i < [expr $nframe - 100]} {incr i} {

    # puts -nonewline [format "      Now is calculating %i_th frames.....\n" $i]
    
    $refsel frame $i
    $refsel update
    
    foreach j [$refsel get index] {

      set sel_tp [atomselect top "index $j"]
      set sel_water [atomselect top "$calatm && pbwithin $rmax of index $j"]
      $sel_water frame $i
      $sel_water update
      
      foreach k [$sel_water get index] {

        set sel_OT [atomselect top "index $k"]
               
        for {set n $i} {$n <= [expr $i + 50]} {incr n} {

          $sel_OT frame $n
          $sel_tp frame $n

          $sel_OT update
          $sel_tp update

          set center_OT [measure center $sel_OT]
          set center_tp [measure center $sel_tp]

          if {$n == $i} {

            set dis 0
            set dis_initial [veclength [vecsub $center_OT $center_tp]]
            
          }
                    
          if {$n != $i} {

            set distmp [expr ([veclength [vecsub $center_OT $center_tp]] - $dis_initial)**2]        
            append dis " "
            append dis $distmp

          }   

        }

        $sel_OT delete
        puts $fileid "$dis" 

      }

      $sel_tp delete
      $sel_water delete

    }
  
  }

  $refsel delete
  close $fileid
  
}

proc integrate_MSD {n} {

  array set data {}
  array unset data *

  set outfile [open "Data_MSD.dat" w]

  set num_cols 0

  for {set i 0} {$i < $n} {incr i} {

    set fd [open Data_MSD_$i\.dat r]

    while {[gets $fd line] >= 0} {

      set values [regexp -all -inline {\S+} $line]
      set current_cols [llength $values]

      if {$num_cols == 0} {

        set num_cols $current_cols

        for {set col 0} {$col < $num_cols} {incr col} {

          set data($col) [list]

        }
        
      }

      if {[llength $values] != $num_cols} {continue}

      for {set col 0} {$col < $num_cols} {incr col} {

        set val [lindex $values $col]
        lappend data($col) $val

      }

    }

    close $fd
        
  }

  set averages [list]

  for {set col 0} {$col < $num_cols} {incr col} {

    set sum 0.0

    set count [llength $data($col)]

    if {$col ==0} {puts $outfile "$count"}

    foreach val $data($col) {

      set sum [expr {$sum + $val}]

    }

    set avg [expr {$count > 0 ? $sum / $count : 0.0}]

    lappend averages $avg

  }

  puts $outfile [join $averages " "]
  close $outfile
  
}

proc measure_MRT {outputname refatm {calatm "type OT"} {rmax "3.5"}} {

  set fileid [open $outputname\.dat w]
  set nframe [molinfo top get numframes]
  set refsel [atomselect top "$refatm"]

  for {set i 1} {$i < [expr $nframe - 100]} {incr i} {

    # puts -nonewline [format "      Now is calculating %i_th frames.....\n" $i]
    
    $refsel frame $i
    $refsel update
    
    foreach j [$refsel get index] {

      set sel_tp [atomselect top "index $j"]
      set sel_water [atomselect top "$calatm && pbwithin $rmax of index $j"]
      $sel_water frame $i
      $sel_water update
      
      foreach k [$sel_water get index] {

        
        set sel_OT [atomselect top "index $k"]
        
        $sel_tp frame [expr $i - 1]
        $sel_OT frame [expr $i - 1]
        $sel_tp update
        $sel_OT update
        set n $i

        unset -nocomplain mark
        set mark  [lrepeat 51 1]
        set logo  0
        set found 0

        if {[veclength [vecsub [measure center $sel_OT] [measure center $sel_tp]]] < $rmax} {

          $sel_OT delete
          continue

        }
        
        for { } {$n <= [expr $i + 50]} {incr n} {

          $sel_OT frame $n
          $sel_tp frame $n
          $sel_OT update
          $sel_tp update
          set distmp [veclength [vecsub [measure center $sel_OT] [measure center $sel_tp]]]

          if {$distmp > $rmax || $found} {
            lset mark $logo 0
            set found 1 

          }

          incr logo

        }

        $sel_OT delete
        puts $fileid [join $mark " "]
         
      }

      $sel_tp delete
      $sel_water delete

    }
  
  }

  $refsel delete
  close $fileid
  
}

proc integrate_MRT {n} {

  array set column_sums {}
  set total_lines 0
  set num_cols 0

  set outfile [open "Data_MRT.dat" w]

  for {set i 0} {$i < $n} {incr i} {

    set fd [open Data_MRT_$i\.dat r]

    while {[gets $fd line] != -1} {

      set values [regexp -all -inline {\S+} $line]

      if {$num_cols == 0} {

        set num_cols [llength $values]
        
        for {set col 0} {$col < $num_cols} {incr col} {

          set column_sums($col) 0.0

        }

      }

      if {[llength $values] != $num_cols} {continue}
      
      for {set col 0} {$col < $num_cols} {incr col} {

        set val [lindex $values $col]
        set column_sums($col) [expr {$column_sums($col) + $val}]

      }
      
      incr total_lines

    }
    
    close $fd
        
  }

  puts $outfile $total_lines

  set averages {}

  for {set col 0} {$col < $num_cols} {incr col} {

    set avg [expr {$total_lines > 0 ? $column_sums($col) / $total_lines : 0.0}]
    lappend averages [format "%.4f" $avg]

  }
  
  puts $outfile [join $averages " "]
  close $outfile
  
}