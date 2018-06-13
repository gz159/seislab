use POSIX;
use Parallel::ForkManager;


# maximum P wave velocity perturbation allowed
# you need to change this value according to you need


$pwd=qx{pwd};
chomp($pwd);

$pm=new Parallel::ForkManager(6);

@dirnames=qx{ls};
foreach $dn(@dirnames){
    
   my $pid=$pm->start and next;
   
   chomp($dn);
   chdir($pwd);
   
   if(-d $dn){
      chdir($dn);
      qx{cp ../Qc.py .};
      qx{python Qc.py};
      qx{rm Qc.py};
      print "Finish process directory: $dn\n";
   }
   $pm->finish;
}
$pm->wait_all_children;
print "******** Finish proces all data *********\n";