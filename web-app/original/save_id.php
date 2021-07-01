<?php
if(!empty($_POST['data'])){
$data = $_POST['data'];
$fname = $_POST['id'];
//$fname = mktime() . ".txt";//generates random name

$file = fopen("output/" .$fname, 'w');//creates new file
fwrite($file, $data);
fclose($file);
}
?>
