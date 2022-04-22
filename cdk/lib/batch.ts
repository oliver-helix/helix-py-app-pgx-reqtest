import * as batch from "@aws-cdk/aws-batch";
import * as ecr from "@aws-cdk/aws-ecr";
import * as ecs from '@aws-cdk/aws-ecs';
import * as iam from '@aws-cdk/aws-iam';
import * as cdk from '@aws-cdk/core';
import { core } from '@myhelix/cdk-library';
import * as ec2 from '@aws-cdk/aws-ec2';
import * as route53  from '@aws-cdk/aws-route53'
import * as s3  from '@aws-cdk/aws-s3'
import {getResourceArn} from "@aws-cdk/aws-stepfunctions-tasks/lib/resource-arn-suffix";

export interface BatchProps {
    vpc: ec2.IVpc,
    bucket: string,
    ecrRepository: ecr.IRepository,
    imageTag: string,
    r2vVersion: string,
    projectName: string,
    environment: core.NamedEnv,
    route53record: string
}

export class Batch extends cdk.Construct {
    readonly jobQueue: batch.JobQueue;
    readonly jobDefinition: batch.JobDefinition;

    constructor(scope: cdk.Construct, id: string, props: BatchProps) {
        super(scope, id);

        const batchJobRole = new iam.Role(this, 'JobRole', {
            assumedBy: new iam.ServicePrincipal('ecs-tasks.amazonaws.com'),
            managedPolicies: [iam.ManagedPolicy.fromAwsManagedPolicyName('service-role/AmazonECSTaskExecutionRolePolicy')],
        });

        const resultsBucket = s3.Bucket.fromBucketName(this, 'r2v_bucket', props.bucket);

        resultsBucket.grantReadWrite(batchJobRole);

        //below bucket only exists in master account
        if (props.environment.account == core.Environment.MasterProduction.account) {
            const awBucketName = `helix-analysis-workflow-results-${props.environment.name}`;
            const awBucket = s3.Bucket.fromBucketName(this, 'aw_bucket', awBucketName);
            awBucket.grantReadWrite(batchJobRole)
        }

        const dataBucket = s3.Bucket.fromBucketName(this, 'data_bucket', 'helix-data-r2v');
        dataBucket.grantRead(batchJobRole)

        const mystack = cdk.Stack.of(this) as core.Stack;

        const jobDefinition = new batch.JobDefinition(this, 'BatchJob', {
            jobDefinitionName: `r2v-batch-job-${mystack.namedEnv.name}`,
            container: {
                image: ecs.ContainerImage.fromEcrRepository(props.ecrRepository, props.imageTag),
                vcpus: 32,
                memoryLimitMiB: 66560, //65 GiB
                jobRole: batchJobRole,
                command: ['r2v', '--platform', 'host', '--json', 'Ref::input_json_s3_url'],
                volumes: [{
                    name: 'ssd',
                    host: {
                        sourcePath: '/scratch'
                    }},
                    {
                    name: 'efs',
                    host: {
                        sourcePath: '/mnt/efs'
                    }}],
                mountPoints: [{
                    containerPath: '/reference',
                    readOnly: false,
                    sourceVolume: 'efs'
                },
                {
                    containerPath: '/scratch',
                    readOnly: false,
                    sourceVolume: 'ssd'
                }],
                environment: {
                    AWS_ENVIRONMENT_NAME: props.environment.name
                }
            },
            retryAttempts: 2,
            timeout: cdk.Duration.hours(6)
        });

        //This is so that when cdk destroy is run, this job definition will remain
        const resource = jobDefinition.node.findChild('Resource') as cdk.CfnResource;
        resource.applyRemovalPolicy(cdk.RemovalPolicy.RETAIN);

        //Route53 record
        let [domainName, recordName] = Batch.getRoute53RecordInfo(mystack.namedEnv, props.route53record)
        const zone =  route53.HostedZone.fromLookup(this, 'myzone', {domainName});
        const route53Record = new route53.TxtRecord(this, 'TXTRecord', {
            zone,
            recordName,
            values: [jobDefinition.jobDefinitionArn, props.r2vVersion],
        });

        this.jobDefinition = jobDefinition;
    }

    static getRoute53RecordInfo(namedEnv: core.NamedEnv, recordRoot: string): [string, string] {
    // Suffix record name with environment name in master account
    if (namedEnv.account == core.Environment.MasterProduction.account) {
        return ['helix.com', `${recordRoot}.${namedEnv.name}`]
    }
    else {
        return [`${namedEnv.name}.helix.com`, recordRoot]
    }

  }
}
