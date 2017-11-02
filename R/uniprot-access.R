
get_uniprot <- function(focal_taxon){

}


    # args = parser()
    #
    # if args.print_http:
    #     httplib2.debuglevel = 1
    #
    # proteome = 'reference:yes' if args.reference else 'complete:yes'
    # include = 'yes' if args.include_isoforms else 'no'
    #
    # h = httplib2.Http(args.cache)
    #
    # par = {'query': 'ancestor:%s+%s' % (args.top_taxid, proteome),
    #        'format': 'list'}
    # response, content = h.request(makeurl('taxonomy', par))
    #
    # if args.print_http:
    #     prettyprint_http(response)
    #
    # taxids = content.decode().strip().split("\n")
    #
    # for taxid in taxids:
    #     if args.retrieve_sequence or args.retrieve_annotations:
    #         if args.retrieve_annotations:
    #             filename = '%s.tab' % taxid
    #             par = {'query':'organism:%s' % taxid,
    #                    'format':'tab',
    #                    'include':include,
    #                    'columns':','.join(args.retrieve_annotations)}
    #         elif args.retrieve_sequence:
    #             filename = '%s.faa' % taxid
    #             par = {'query':'organism:%s' % taxid,
    #                    'format':'fasta',
    #                    'include':include}
    #         print("Retrieving %s" % filename, file=sys.stderr)
    #         url = makeurl('uniprot', par)
    #         response, content = h.request(url)
    #         if args.print_http:
    #             print("URL:%s" % url)
    #             prettyprint_http(response)
    #         with open(filename, 'w') as f:
    #             print(content.decode(), file=f)
    #     else:
    #         print(taxid)
