function table_to_csv(table) {
    const columns = Array.from(table.columns, ({ field }) => field)
    const nrows = table.source.get_length()
    const lines = [columns.join(',')]

    console.log(lines)
    console.log(nrows)
    console.log(table.view.filters[0])

    for (let i = 0; i < nrows; i++) {
        if (table.view.filters[0].booleans[i]) {
            let row = [];
            for (let j = 0; j < columns.length; j++) {
                const column = columns[j]
                let cell = table.source.data[column][i].toString()
                if (cell.includes(',')) {
                    cell = '"' + cell + '"'
                }
                row.push(cell)
            }
            lines.push(row.join(','))
        }
    }
    return lines.join('\n').concat('\n')
}

const filename = 'data_result.csv'
console.log(filename)
const filetext = table_to_csv(table)
const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}