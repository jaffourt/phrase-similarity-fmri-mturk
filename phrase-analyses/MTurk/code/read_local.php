<?php
    $file = fopen('data.csv', 'r');
    $string = "";
    while (($line = fgetcsv($file)) !== FALSE) {
        $string .= $line[0] . "," . $line[1] . "," . $line[2] . "\n";
    }
    fclose($file);
    echo $string
?>
