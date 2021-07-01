$.ajax({
    type: 'get',
    url: 'read_local.php',
    data: '',
    success: function(data) {
    alert(data);
       }
});
//var data = $.csv.toObjects(csv)

