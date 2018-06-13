use POSIX;
use Parallel::ForkManager;


# maximum P wave velocity perturbation allowed
# you need to change this value according to you need


$pwd=qx{pwd};
chomp($pwd);

#$pm=new Parallel::ForkManager(6);

@dirnames=qx{ls};
foreach $dn(@dirnames){
    
   #my $pid=$pm->start and next;
   
   chomp($dn);
   chdir($pwd);
   
   
   if(-d $dn){
    chdir($dn);
    qx{rm *.psmov};
    qx{rm *.stk};
    @dnf=split('\.',$dn);
    $sacfn=R.".".$dnf[1].".".$dnf[0].".*".".1.5";

    @filenames=qx{ls $sacfn};
    foreach $fn(@filenames){
      chomp($fn);
      $user0=qx{saclhdr -USER0 $fn};
      chomp($user0);
      
      #pythyon calculate the flowing number
      #>>> from obspy.geodetics import kilometer2degrees
      #>>> kilometer2degrees(1)
      $slowness=$user0/0.008993216059187306;
      chomp($slowness);

      # add SAC header value nzsec to 00
      # otherwise the read_rf will compline
      # SAC file header value error
      open(SAC, "|sac\n");
      print SAC "read $fn\n";
      print SAC "ch nzsec 00\n";
      print SAC "wh\n";
      print SAC "write over\n";
      print SAC "quit\n";
      close(SAC);

      $outfn=$fn.".psmov";
      chomp($outfn);

      open(PY, "|python\n");
      print PY "from rf import read_rf\n";
      print PY "s=read_rf(\'$fn\')\n";
      print PY "s[0].stats.inclination = $slowness\n";
      print PY "s[0].stats.onset = s[0].stats.event_time\n";
      print PY "s.moveout(phase='Ps')\n";
      print PY "s.write(\'$outfn\',\'SAC\')\n";
      close(PY);
      
      # change to the original user0 value
      open(SAC, "|sac\n");
      print SAC "read $outfn\n";
      print SAC "ch user0 $user0\n";
      print SAC "wh\n";
      print SAC "write over\n";
      print SAC "quit\n";
      close(SAC);
      
    }
    
    # plot Ps move out correction rf and stacking result 
    open(PY, "|python\n");
    print PY "from rf import read_rf\n";
    print PY "s=read_rf(\'*.psmov\')\n";
    print PY "s.plot_rf(fname=\'sktrf.png\',trace_height=0.30, stack_height=0.5)\n";
    close(PY);
    
    #open(SAC, "|gsac\n");
    #print SAC "read *.psmov\n";
    #print SAC "stack relative norm on\n";
    #print SAC "cd ..\n";
    #print SAC "write\n";
    #print SAC "quit\n";
    #close(SAC);

    
   }
   #$pm->finish;
}
#$pm->wait_all_children;
